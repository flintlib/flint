/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "gr.h"
#include "nfloat.h"
#include "impl.h"
#include "elem.h"

/*
    Square root: for x = A 2^e with A = a / B^n in [1/2, 1) and t = e mod 2,
    sqrt(x) = sqrt(A / 2^t) 2^((e + t) / 2) where sqrt(A / 2^t) is in
    [1/2, 1). Its truncated mantissa floor(sqrt(a B^n / 2^t)) is computed
    exactly by flint_mpn_sqrtrem, which also tells whether the root is
    exact, so all rounding modes are handled exactly.
*/
int
nfloat_sqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    ulong T[2 * NFLOAT_MAX_LIMBS];
    nn_srcptr a;
    slong n, e, exp;
    ulong t;
    int inexact;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_NEG_INF(x))
            return nfloat_nan(res, ctx);
        else
            return nfloat_set(res, x, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    t = e & 1;
    a = NFLOAT_D(x);

    if (n == 1)
    {
#if FLINT_BITS == 64
        /* the exact floor of sqrt(a 2^(FLINT_BITS - t)), in registers */
        ulong s1;

        inexact = _nfloat_isqrt_2(&s1, a[0] >> t, (a[0] << (FLINT_BITS - 1)) & (-t));
        exp = (e + (slong) t) / 2;
        NFLOAT_D(res)[0] = s1;
        if (inexact && nfloat_should_round_up(0, ctx))
            NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);
        NFLOAT_EXP(res) = exp;
        NFLOAT_SGNBIT(res) = 0;
        return GR_SUCCESS;
#else
        T[1] = a[0] >> t;
        T[0] = (a[0] << (FLINT_BITS - 1)) & (-t);
#endif
    }
#if FLINT_BITS == 64
    else if (n == 2)
    {
        /* sqrt(A 2^-t) at three limbs within NFLOAT_SMALL_SQRT_ERR_3
           units: the truncation to two limbs is certified, and inexact,
           when the guard limb is not within that many units of 0 or B */
        ulong s2, s1, s0;

        _nfloat_small_sqrt_3(&s2, &s1, &s0, a[1], a[0], 0, t);
        if (s0 > NFLOAT_SMALL_SQRT_ERR_3 && s0 < UWORD_MAX - NFLOAT_SMALL_SQRT_ERR_3)
        {
            exp = (e + (slong) t) / 2;
            NFLOAT_D(res)[1] = s2;
            NFLOAT_D(res)[0] = s1;
            if (nfloat_should_round_up(0, ctx))
                NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);
            NFLOAT_EXP(res) = exp;
            NFLOAT_SGNBIT(res) = 0;
            return GR_SUCCESS;
        }
        _nfloat_zero_limbs(T, n - 1);
        if (t)
            T[n - 1] = mpn_rshift(T + n, a, n, 1);
        else
        {
            T[n - 1] = 0;
            _nfloat_copy_limbs(T + n, a, n);
        }
    }
#endif
    else
    {
        _nfloat_zero_limbs(T, n - 1);
        if (t)
        {
            T[n - 1] = mpn_rshift(T + n, a, n, 1);
        }
        else
        {
            T[n - 1] = 0;
            _nfloat_copy_limbs(T + n, a, n);
        }
    }

    inexact = (flint_mpn_sqrtrem(NFLOAT_D(res), NULL, T, 2 * n) != 0);
    exp = (e + (slong) t) / 2;

    if (inexact && nfloat_should_round_up(0, ctx))
        NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);

    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}

/*
    Reciprocal square root: with x = A 2^e, t = e mod 2 and
    M = (4A)^(-1/2) in (1/2, 2^(-1/2)] (t = 0) or M = (2A)^(-1/2) in
    (2^(-1/2), 1) (t = 1, A != 1/2), we have
    1/sqrt(x) = M 2^((1 - t) - (e - t) / 2).

    Exactly, floor(M B^n) = floor(sqrt(floor(z))) with z = B^(3n) / (4a)
    resp. B^(3n) / (2a), and M B^n is an integer iff both the division
    and the square root are exact.
*/

/* (s, n) = floor(M B^n); returns inexactness. s may alias a. */
static int
_nfloat_rsqrt_mpn_exact(nn_ptr s, nn_srcptr a, slong n, int t)
{
    ulong N[3 * NFLOAT_MAX_LIMBS];
    ulong Q[2 * NFLOAT_MAX_LIMBS + 1];
    ulong R[NFLOAT_MAX_LIMBS];
    int inexact;

    _nfloat_zero_limbs(N, 3 * n - 1);
    N[3 * n - 1] = UWORD(1) << (FLINT_BITS - 2 + t);
    flint_mpn_tdiv_qr(Q, R, N, 3 * n, a, n);
    FLINT_ASSERT(Q[2 * n] == 0);
    inexact = !flint_mpn_zero_p(R, n);
    inexact |= (flint_mpn_sqrtrem(s, NULL, Q, 2 * n) != 0);
    return inexact;
}


/*
    Newton iteration for n >= 2 limbs, computing (Y, n + 1) with
    |Y - M B^(n+1)| <= E (E = 196 if n = 3, E = 1.03 otherwise), so that
    the n-limb truncation is certified (and inexact) when the guard limb
    D = Y[0] satisfies E < D < B - E.

    Base: two limbs from the top two limbs (a1, a0) of a, as in the
    one-limb kernel above but with the residual 2^(254 + t) - a Y0^2
    computed from both limbs; the truncation of a changes M by at most
    one unit of B^-2, so |Y - M B^2| < 8.

    Step (Y, m) -> (Z, p) with m < p <= 2m, given |y - M| <= e_m where
    y = Y / B^m: with L = p + 2, h is the mulhigh of the top L limbs of
    a and of S = Y^2 (both lower bounds of A resp. s = y^2, the product a
    lower bound with error below 2 units of B^-L), so that
    0 <= A s - h < 4 B^-L. The residual 1 - c y^2 is taken as
    1 - 2^(2-t) h, an overestimate by less than 16 B^-L, truncated to
    p + 1 fraction limbs (error < B^-(p+1)). The correction y eps / 2 is
    truncated to p fraction limbs (error < B^-p). Exact Newton steps
    have error at most 1.5 eta^2 M (1 + eta) <= 3.04 e_m^2 (M > 1/2).
    Hence e_p <= 3.04 e_m^2 + 8 B^-(p+2) + 0.5 B^-(p+1) + B^-p, i.e. with
    e_m = C B^-m, e_p <= (3.04 C^2 B^(p - 2m) + 1.01) B^-p. The first step
    (C = 8) may go up to p = 2m = 4, giving e_4 < 196 B^-4; all further
    steps have p <= 2m - 1, giving e_p < 1.03 B^-p. Results that would
    reach B^p are clamped to B^p - 1 (M B^p < B^p, so the bound is
    preserved).
*/

#if FLINT_BITS == 64

/* (Y, 2) from the top two limbs of a (elem.h) */
static void
_nfloat_rsqrt_newton_base(nn_ptr Y, ulong a1, ulong a0, ulong t)
{
    _nfloat_small_rsqrt_2(Y + 1, Y, a1, a0, t);
}

#else

/* (Y, 2) = floor(M' B^2) for the truncation a' of a to two limbs, computed
   exactly (clamped to B^2 - 1 when M' = 1); |Y - M B^2| < 2 */
static void
_nfloat_rsqrt_newton_base(nn_ptr Y, ulong a1, ulong a0, ulong t)
{
    ulong a2[2];

    if (t && a1 == (UWORD(1) << (FLINT_BITS - 1)) && a0 == 0)
    {
        Y[0] = Y[1] = UWORD_MAX;
        return;
    }

    a2[0] = a0;
    a2[1] = a1;
    _nfloat_rsqrt_mpn_exact(Y, a2, 2, t);
}

#endif

/* (Z, p) from (Y, m), m < p <= 2m */
static void
_nfloat_rsqrt_newton_step(nn_ptr Z, slong p, nn_srcptr Y, slong m, nn_srcptr a, slong n, ulong t)
{
    ulong S[2 * NFLOAT_MAX_LIMBS + 4];
    ulong ap[NFLOAT_MAX_LIMBS + 4];
    ulong Sp[NFLOAT_MAX_LIMBS + 4];
    ulong H[NFLOAT_MAX_LIMBS + 4];
    ulong G[NFLOAT_MAX_LIMBS + 5];
    ulong T[NFLOAT_MAX_LIMBS + 6];
    nn_srcptr ah, Sh;
    nn_ptr E1;
    slong L = p + 2, en, dn;
    int neg;

    FLINT_ASSERT(m >= 2 && m < p && p <= 2 * m);

    flint_mpn_sqr(S, Y, m);

    if (n >= L)
        ah = a + n - L;
    else
    {
        _nfloat_zero_limbs(ap, L - n);
        _nfloat_copy_limbs(ap + L - n, a, n);
        ah = ap;
    }

    if (2 * m >= L)
        Sh = S + 2 * m - L;
    else
    {
        _nfloat_zero_limbs(Sp, L - 2 * m);
        _nfloat_copy_limbs(Sp + L - 2 * m, S, 2 * m);
        Sh = Sp;
    }

    flint_mpn_mulhigh_n(H, ah, Sh, L);

    /* G = 2^(2 - t) h B^L ~ B^L; |E| = |B^L - G| < B^(L - m + 1) */
    G[L] = mpn_lshift(G, H, L, 2 - t);
    neg = (G[L] != 0);
    FLINT_ASSERT(G[L] <= 1);
    if (!neg)
        mpn_neg(G, G, L);
    FLINT_ASSERT(flint_mpn_zero_p(G + L - m + 1, m - 1));

    /* E1 = floor(|E| / B), p + 2 - m limbs */
    E1 = G + 1;
    en = L - m;
    while (en > 0 && E1[en - 1] == 0)
        en--;

    _nfloat_zero_limbs(Z, p - m);
    _nfloat_copy_limbs(Z + p - m, Y, m);

    if (en == 0)
        return;

    /* T = Y E1 (at most p + 2 limbs); D = T / (2 B^(m + 1)) */
    if (m >= en)
        flint_mpn_mul(T, Y, m, E1, en);
    else
        flint_mpn_mul(T, E1, en, Y, m);
    _nfloat_zero_limbs(T + m + en, p + 2 - m - en);

    dn = p + 1 - m;
    mpn_rshift(T + m + 1, T + m + 1, dn, 1);

    if (neg)
    {
        mpn_sub(Z, Z, p, T + m + 1, dn);
    }
    else
    {
        if (mpn_add(Z, Z, p, T + m + 1, dn))
            flint_mpn_store(Z, p, UWORD_MAX);
    }
}

/* (Y, n + 1) with |Y - M B^(n + 1)| <= 196 (n = 3) resp. 1.03, n >= 2 */
static void
_nfloat_rsqrt_newton(nn_ptr Y, nn_srcptr a, slong n, ulong t)
{
    ulong W[2][NFLOAT_MAX_LIMBS + 2];
    slong prec[FLINT_BITS];
    slong i, k;
    int cur;

    k = 0;
    prec[0] = n + 1;
    while (prec[k] > 4)
    {
        prec[k + 1] = (prec[k] + 2) / 2;
        k++;
    }
    prec[k + 1] = 2;
    k++;

    cur = 0;
    _nfloat_rsqrt_newton_base(W[cur], a[n - 1], a[n - 2], t);

    for (i = k - 1; i >= 0; i--)
    {
        _nfloat_rsqrt_newton_step((i == 0) ? Y : W[cur ^ 1], prec[i], W[cur], prec[i + 1], a, n, t);
        cur ^= 1;
    }
}

int
nfloat_rsqrt(nfloat_ptr res, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n, e, exp;
    ulong t;
    nn_srcptr a;
    int inexact;

    if (NFLOAT_IS_SPECIAL(x))
    {
        if (NFLOAT_IS_ZERO(x))
            return nfloat_pos_inf(res, ctx);
        else if (NFLOAT_IS_POS_INF(x))
            return nfloat_zero(res, ctx);
        else
            return nfloat_nan(res, ctx);
    }

    if (NFLOAT_SGNBIT(x))
        return nfloat_nan(res, ctx);

    n = NFLOAT_CTX_NLIMBS(ctx);
    e = NFLOAT_EXP(x);
    t = e & 1;
    a = NFLOAT_D(x);

    /* x = 2^(e - 1) with e - 1 even */
    if (t && a[n - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(a, n - 1))
    {
        if (res != x)
            _nfloat_copy_limbs(NFLOAT_D(res), a, n);
        NFLOAT_EXP(res) = (1 - e) / 2 + 1;
        NFLOAT_SGNBIT(res) = 0;
        return GR_SUCCESS;
    }

    exp = (1 - (slong) t) - (e - (slong) t) / 2;

    if (n == 1)
    {
#if FLINT_BITS == 64
        ulong m;
        int directed = NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx);

        if (_nfloat_rsqrt_1_approx(&m, a[0], t, directed))
        {
            NFLOAT_D(res)[0] = m;
            inexact = 1;
        }
        else
#endif
        {
            inexact = _nfloat_rsqrt_mpn_exact(NFLOAT_D(res), a, 1, t);
        }
    }
#if FLINT_BITS == 64
    else if (n == 2)
    {
        /* three limbs in registers within NFLOAT_SMALL_RSQRT_ERR_3 units */
        ulong z2, z1, z0;

        _nfloat_small_rsqrt_3(&z2, &z1, &z0, a[1], a[0], t);

        if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx) ||
            (z0 > NFLOAT_SMALL_RSQRT_ERR_3 && z0 < UWORD_MAX - NFLOAT_SMALL_RSQRT_ERR_3))
        {
            NFLOAT_D(res)[1] = z2;
            NFLOAT_D(res)[0] = z1;
            inexact = NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx);
        }
        else
        {
            inexact = _nfloat_rsqrt_mpn_exact(NFLOAT_D(res), a, n, t);
        }
    }
#endif
    else
    {
        ulong Y[NFLOAT_MAX_LIMBS + 1];

        _nfloat_rsqrt_newton(Y, a, n, t);

        if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        {
            _nfloat_copy_limbs(NFLOAT_D(res), Y + 1, n);
            inexact = 0;
        }
        else if ((n == 3) ? (Y[0] >= 197 && Y[0] <= UWORD_MAX - 197)
                          : (Y[0] >= 2 && Y[0] <= UWORD_MAX - 2))
        {
            _nfloat_copy_limbs(NFLOAT_D(res), Y + 1, n);
            inexact = 1;
        }
        else
        {
            inexact = _nfloat_rsqrt_mpn_exact(NFLOAT_D(res), a, n, t);
        }
    }

    if (inexact && nfloat_should_round_up(0, ctx))
        NFLOAT_MANT_INCREMENT(NFLOAT_D(res), n, exp);

    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = 0;
    return GR_SUCCESS;
}
