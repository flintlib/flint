/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "fp_small.h"

/*
    Floating-point square roots and reciprocal square roots with a guard
    limb (see mpn_extras.h). With A = {a, an} / B^an in [1/2, 1) and t the
    parity of the exponent of x = A 2^e, the quantities are

        S = sqrt(A / 2^t) in [1/2, 1),     sqrt(x) = S 2^((e + t) / 2),
        M = (2^(2-t) A)^(-1/2) in (1/2, 1], 1/sqrt(x) = M 2^((1 - t) - (e - t) / 2),

    and the n + 1 limbs approximate S B^(n+1) resp. M B^(n+1).
*/

/* (X, xn = 2 (n + 1)) = floor(A B^xn / 2^t) with A truncated to xn limbs;
   returns the dropped part (limbs below and the bit shifted out) as a
   nonzero flag, and with cmp != NULL, the sign of f - 1/4 for the
   dropped fraction f = A B^xn / 2^t - X in [0, 1) */
FLINT_FORCE_INLINE ulong
_sqrt_operand(nn_ptr X, int * cmp, nn_srcptr a, slong an, int t, slong n)
{
    slong xn = 2 * (n + 1);
    ulong dropped = 0;
    int c = -1;

    if (an >= xn)
    {
        slong k = an - xn;
        ulong bit = 0;

        if (t)
            bit = mpn_rshift(X, a + k, xn, 1) >> (FLINT_BITS - 1);
        else
            _flint_mpn_copy_small(X, a + k, xn);

        dropped = bit | !flint_mpn_zero_p(a, k);

        if (cmp != NULL && dropped)
        {
            /* f = bit / 2 + L / (2^t B^k) for the k dropped limbs L */
            if (bit)
                c = 1;
            else if (k > 0)
            {
                ulong thresh = UWORD(1) << (FLINT_BITS - 1 - (t == 0));
                if (a[k - 1] != thresh)
                    c = (a[k - 1] > thresh) ? 1 : -1;
                else
                    c = !flint_mpn_zero_p(a, k - 1);
            }
        }
    }
    else
    {
        /* the bit shifted out lands in the padding */
        slong off = xn - an;
        _flint_mpn_zero_small(X, off);
        if (t)
            X[off - 1] = mpn_rshift(X + off, a, an, 1);
        else
            _flint_mpn_copy_small(X + off, a, an);
    }

    if (cmp != NULL)
        *cmp = c;
    return dropped;
}

/* operands of up to this many limbs are held on the stack */
#define SQRTHIGH_STACK_LIMBS 64

/* for the Newton routines of mp_real: the fraction A / 2^t in [1/4, 1)
   truncated to m limbs (m >= 2), written to c (m limbs) */
static void
_sqrt_operand_fraction(nn_ptr c, nn_srcptr a, slong an, int t, slong m)
{
    if (an >= m)
        flint_mpn_copyi(c, a + an - m, m);
    else
    {
        flint_mpn_zero(c, m - an);
        flint_mpn_copyi(c + m - an, a, an);
    }
    if (t)
        mpn_rshift(c, c, m, 1);
}

/* the square root of an integer of up to this many limbs is taken by
   flint_mpn_sqrtrem on the zero-extended operand; above, by the Newton
   iteration of mp_real, which needs no remainder (as mp_real_sqrt) */
#define SQRTHIGH_NEWTON_CUTOFF (FLINT_MPN_SQRTREM_NEWTON_CUTOFF / 4)

void
flint_mpn_sqrthigh(nn_ptr s, nn_srcptr a, slong an, int t, slong n)
{
    TMP_INIT;

    FLINT_ASSERT(n >= 1 && an >= 1);
    FLINT_ASSERT(a[an - 1] >> (FLINT_BITS - 1));

#if FLINT_BITS == 64
    if (n == 1)
    {
        /* two limbs from the top two limbs of A (truncation: one unit) */
        ulong a1 = a[an - 1], a0 = (an >= 2) ? a[an - 2] : 0, s1, s0;
        _flint_mpn_small_sqrt_2(&s1, &s0, a1, a0, t);
        s[1] = s1;
        s[0] = s0;
        return;
    }
    if (n == 2)
    {
        ulong a2 = a[an - 1], a1 = (an >= 2) ? a[an - 2] : 0, a0 = (an >= 3) ? a[an - 3] : 0, s2, s1, s0;
        _flint_mpn_small_sqrt_3(&s2, &s1, &s0, a2, a1, a0, t);
        s[2] = s2;
        s[1] = s1;
        s[0] = s0;
        return;
    }
#endif

    TMP_START;

    if (n + 1 < SQRTHIGH_NEWTON_CUTOFF)
    {
        /* floor(sqrt(X)) for X = floor(A B^(2n+2) / 2^t): the root of the
           truncated X is below the exact value by less than
           1 + B^(-2n-2) / (2 S) < 1 + 2 / B */
        slong xn = 2 * (n + 1);
        ulong stack[SQRTHIGH_STACK_LIMBS];
        nn_ptr X = (xn <= SQRTHIGH_STACK_LIMBS) ? stack : TMP_ALLOC(xn * sizeof(ulong));

        _sqrt_operand(X, NULL, a, an, t, n);
        flint_mpn_sqrtrem(s, NULL, X, xn);
    }
    else
    {
        /* sqrt(c) for c = A / 2^t in [1/4, 1) at n + 2 fraction limbs
           (error 4 B^-(n+2) / sqrt(c) <= 8 B^-(n+2)), the lowest limb
           dropped */
        slong m = n + 4, nd = n + 2;
        nn_ptr c = TMP_ALLOC((m + nd + 2) * sizeof(ulong)), q = c + m;

        _sqrt_operand_fraction(c, a, an, t, m);
        _mp_real_sqrt_newton(q, c, m, nd);
        FLINT_ASSERT(q[nd] == 0 && q[nd + 1] == 0);
        flint_mpn_copyi(s, q + 1, n + 1);
    }

    /* the top bit is set unless the approximation dipped below 1/2 */
    if (!(s[n] >> (FLINT_BITS - 1)))
    {
        flint_mpn_zero(s, n);
        s[n] = UWORD(1) << (FLINT_BITS - 1);
    }

    TMP_END;
}

/*
    The exact floor, with the inexact flag alone (the remainder of the
    integer root not formed, which is cheaper for long inputs), with
    the remainder, or with the position of the root relative to the
    next half unit: with X = floor(A B^xn / 2^t) the radicand of the
    n + 1 limbs (xn = 2 (n + 1)), f in [0, 1) the dropped part of A and
    R = X - F^2 the remainder of the integer root F, the root exceeds
    F + 1/2 iff X + f >= F^2 + F + 1/4, i.e. R > F, or R = F and
    f >= 1/4 (a tie, f = 1/4, needs dropped limbs).
*/
int
flint_mpn_sqrthigh_exact(nn_ptr s, nn_srcptr a, slong an, int t, slong n)
{
    slong xn = 2 * (n + 1);
    ulong stack[SQRTHIGH_STACK_LIMBS], dropped;
    nn_ptr X;
    int inexact;

    /* n = 0 (a single limb) is allowed here: nfloat uses it on 32-bit
       limbs, where there is no register root */
    FLINT_ASSERT(n >= 0 && an >= 1);

    X = (xn <= SQRTHIGH_STACK_LIMBS) ? stack : flint_malloc(xn * sizeof(ulong));
    dropped = _sqrt_operand(X, NULL, a, an, t, n);
    inexact = (flint_mpn_sqrtrem(s, NULL, X, xn) != 0) || (dropped != 0);
    if (X != stack)
        flint_free(X);
    return inexact;
}

/* the root with the remainder (R, n + 2), zero-padded; *cmp the sign of
   f - 1/4 and *dropped whether f != 0 */
static void
_sqrthigh_exact_rem(nn_ptr s, nn_ptr R, int * cmp, ulong * dropped, nn_srcptr a, slong an, int t, slong n)
{
    slong xn = 2 * (n + 1), rn, need;
    ulong stack[SQRTHIGH_STACK_LIMBS];
    nn_ptr X, W;

    FLINT_ASSERT(n >= 0 && an >= 1);

    /* GMP's mpn_sqrtrem (used on 32-bit machines) needs xn limbs of
       remainder space, taken after X and copied out */
    need = FLINT_MPN_SQRTREM_USE_GMP(xn) ? 2 * xn : xn;
    X = (need <= SQRTHIGH_STACK_LIMBS) ? stack : flint_malloc(need * sizeof(ulong));
    W = (need == xn) ? R : X + xn;
    /* (_sqrt_operand writes all xn limbs; zeroed first since gcc 15
       cannot see that through its branches) */
    flint_mpn_zero(X, xn);
    *dropped = _sqrt_operand(X, cmp, a, an, t, n);
    rn = flint_mpn_sqrtrem(s, W, X, xn);
    if (W != R)
        flint_mpn_copyi(R, W, rn);
    flint_mpn_zero(R + rn, n + 2 - rn);
    if (X != stack)
        flint_free(X);
}

int
flint_mpn_sqrthigh_exact_rem(nn_ptr s, nn_ptr r, nn_srcptr a, slong an, int t, slong n)
{
    ulong dropped;
    int cmp;

    _sqrthigh_exact_rem(s, r, &cmp, &dropped, a, an, t, n);
    return (dropped != 0) || !flint_mpn_zero_p(r, n + 2);
}

int
flint_mpn_sqrthigh_exact_half(nn_ptr s, nn_srcptr a, slong an, int t, slong n)
{
    ulong stack[SQRTHIGH_STACK_LIMBS / 2 + 2], dropped;
    nn_ptr R;
    int cmp, c;

    R = (n + 2 <= SQRTHIGH_STACK_LIMBS / 2 + 2) ? stack : flint_malloc((n + 2) * sizeof(ulong));
    _sqrthigh_exact_rem(s, R, &cmp, &dropped, a, an, t, n);

    /* R against F (n + 1 limbs, top bit set) */
    if (R[n + 1] != 0)
        c = 1;
    else
        c = mpn_cmp(R, s, n + 1);

    if (c > 0)
        c = 3;
    else if (c == 0)
        c = (cmp > 0) ? 3 : (cmp == 0) ? 2 : 1;
    else
        c = (dropped != 0) || !flint_mpn_zero_p(R, n + 1);

    if (R != stack)
        flint_free(R);
    return c;
}

/*
    Reciprocal square root by Newton iteration, computing (Y, N) with
    |Y - M B^N| <= 196 (N = 4) resp. 1.03 (all other N >= 3), 8 (N = 2).

    Base: two limbs from the top two limbs (a1, a0) of a, the double
    precision 1/sqrt refined at 128 bits with the residual of both limbs
    and the third-order term (fp_small.h); the truncation of a changes M
    by at most one unit of B^-2, so |Y - M B^2| < 8.

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

/* (Y, 2) from the top two limbs of a (fp_small.h) */
static void
_rsqrt_newton_base(nn_ptr Y, ulong a1, ulong a0, ulong t)
{
    _flint_mpn_small_rsqrt_2(Y + 1, Y, a1, a0, t);
}

#else

/* (Y, 2) = floor(M' B^2) for the truncation a' of a to two limbs, computed
   exactly (clamped to B^2 - 1 when M' = 1); |Y - M B^2| < 2 */
static void
_rsqrt_newton_base(nn_ptr Y, ulong a1, ulong a0, ulong t)
{
    ulong a2[2];

    if (t && a1 == (UWORD(1) << (FLINT_BITS - 1)) && a0 == 0)
    {
        Y[0] = Y[1] = UWORD_MAX;
        return;
    }

    a2[0] = a0;
    a2[1] = a1;
    flint_mpn_rsqrthigh_exact(Y, a2, 2, t, 1);
}

#endif

/* copies in the Newton step: a few limbs in registers, longer ones by
   memcpy / memset (the compiler does not reliably make that choice) */
FLINT_FORCE_INLINE void
_rsqrt_copy(nn_ptr d, nn_srcptr s, slong n)
{
    if (n <= 8)
        _flint_mpn_copy_small(d, s, n);
    else
        memcpy(d, s, n * sizeof(ulong));
}

FLINT_FORCE_INLINE void
_rsqrt_zero(nn_ptr d, slong n)
{
    if (n <= 8)
        _flint_mpn_zero_small(d, n);
    else
        memset(d, 0, n * sizeof(ulong));
}

/* (Z, p) from (Y, m), m < p <= 2m; scratch of 8 p + 24 limbs */
static void
_rsqrt_newton_step(nn_ptr Z, slong p, nn_srcptr Y, slong m, nn_srcptr a, slong an, ulong t, nn_ptr scratch)
{
    nn_ptr S = scratch, ap = S + 2 * m, Sp = ap + p + 4, H = Sp + p + 4, G = H + p + 4, T = G + p + 5;
    nn_srcptr ah, Sh;
    nn_ptr E1;
    slong L = p + 2, en, dn;
    int neg;

    FLINT_ASSERT(m >= 2 && m < p && p <= 2 * m);

    flint_mpn_sqr(S, Y, m);

    if (an >= L)
        ah = a + an - L;
    else
    {
        _rsqrt_zero(ap, L - an);
        _rsqrt_copy(ap + L - an, a, an);
        ah = ap;
    }

    if (2 * m >= L)
        Sh = S + 2 * m - L;
    else
    {
        _rsqrt_zero(Sp, L - 2 * m);
        _rsqrt_copy(Sp + L - 2 * m, S, 2 * m);
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

    _rsqrt_zero(Z, p - m);
    _rsqrt_copy(Z + p - m, Y, m);

    if (en == 0)
        return;

    /* T = Y E1 (at most p + 2 limbs); D = T / (2 B^(m + 1)) */
    if (m >= en)
        flint_mpn_mul(T, Y, m, E1, en);
    else
        flint_mpn_mul(T, E1, en, Y, m);
    _rsqrt_zero(T + m + en, p + 2 - m - en);

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

/* above this many result limbs, the Newton iteration of mp_real (built
   on the short products of the FFT range) replaces the one below
   (measured crossover: about 110 limbs on Skylake, 96 on Zen 3) */
#define RSQRTHIGH_MP_REAL_CUTOFF FLINT_MPN_RSQRTHIGH_MP_REAL_CUTOFF

/* scratch of the iteration for N result limbs; held on the stack up to
   RSQRT_NEWTON_STACK_N limbs (about 5 KB) */
#define RSQRT_NEWTON_SCRATCH(N) (2 * ((N) + 1) + 8 * (N) + 24)
#define RSQRT_NEWTON_STACK_N 64
#define RSQRT_NEWTON_STACK_LIMBS RSQRT_NEWTON_SCRATCH(RSQRT_NEWTON_STACK_N)

/* (Y, N), N >= 2 (exported for tuning) */
void
_flint_mpn_rsqrthigh_newton(nn_ptr Y, slong N, nn_srcptr a, slong an, ulong t)
{
    slong prec[FLINT_BITS], i, k;
    ulong stack[RSQRT_NEWTON_STACK_LIMBS];
    nn_ptr W[2], scratch;
    int cur;

    if (N == 2)
    {
        /* the starting value is the result (FLINT_BITS == 32: not
           reached from flint_mpn_rsqrthigh on 64-bit limbs) */
        _rsqrt_newton_base(Y, a[an - 1], (an >= 2) ? a[an - 2] : 0, t);
        return;
    }

    k = 0;
    prec[0] = N;
    while (prec[k] > 4)
    {
        prec[k + 1] = (prec[k] + 2) / 2;
        k++;
    }
    prec[k + 1] = 2;
    k++;

    W[0] = (RSQRT_NEWTON_SCRATCH(N) <= RSQRT_NEWTON_STACK_LIMBS) ? stack
        : flint_malloc(RSQRT_NEWTON_SCRATCH(N) * sizeof(ulong));
    W[1] = W[0] + N + 1;
    scratch = W[1] + N + 1;

    cur = 0;
    _rsqrt_newton_base(W[cur], a[an - 1], (an >= 2) ? a[an - 2] : 0, t);

    for (i = k - 1; i >= 0; i--)
    {
        _rsqrt_newton_step((i == 0) ? Y : W[cur ^ 1], prec[i], W[cur], prec[i + 1], a, an, t, scratch);
        cur ^= 1;
    }

    if (W[0] != stack)
        flint_free(W[0]);
}

void
flint_mpn_rsqrthigh(nn_ptr y, nn_srcptr a, slong an, int t, slong n)
{
    FLINT_ASSERT(n >= 1 && an >= 1);
    FLINT_ASSERT(a[an - 1] >> (FLINT_BITS - 1));

#if FLINT_BITS == 64
    if (n == 1)
    {
        /* |Y - M B^2| < 8, plus the truncation of A to two limbs */
        ulong a1 = a[an - 1], a0 = (an >= 2) ? a[an - 2] : 0;
        _flint_mpn_small_rsqrt_2(y + 1, y, a1, a0, t);
        return;
    }
#endif

    if (n + 1 <= RSQRTHIGH_MP_REAL_CUTOFF)
    {
        _flint_mpn_rsqrthigh_newton(y, n + 1, a, an, t);
    }
    else
    {
        /* 1/sqrt(c) for c = A / 2^t in [1/4, 1) at n + 2 fraction limbs
           (error 4 B^-(n+2) / sqrt(c) <= 8 B^-(n+2)) is 2 M: halved, the
           lowest limb dropped */
        slong m = n + 4, nd = n + 2;
        nn_ptr c, q;
        TMP_INIT;

        TMP_START;
        c = TMP_ALLOC((m + nd + 2) * sizeof(ulong));
        q = c + m;
        _sqrt_operand_fraction(c, a, an, t, m);
        _mp_real_rsqrt_newton(q, c, m, nd);
        mpn_rshift(q, q, nd + 2, 1);
        if (q[nd + 1] != 0 || q[nd] != 0)
            flint_mpn_store(y, n + 1, UWORD_MAX);   /* M = 1 */
        else
            flint_mpn_copyi(y, q + 1, n + 1);
        TMP_END;
    }

    if (!(y[n] >> (FLINT_BITS - 1)))
    {
        flint_mpn_zero(y, n);
        y[n] = UWORD(1) << (FLINT_BITS - 1);
    }
}

/*
    floor(M B^(n+1)) = floor(sqrt(floor(Z))) with Z = B^(2(n+1) + an) /
    (2^(2-t) a), M B^(n+1) being an integer iff both the division and the
    square root are exact. Requires M < 1, i.e. not (t = 1 and A = 1/2).
*/
int
flint_mpn_rsqrthigh_exact(nn_ptr y, nn_srcptr a, slong an, int t, slong n)
{
    slong xn = 2 * (n + 1), nn = xn + an + 1, dn = an + 1;
    nn_ptr N, D, Q, R;
    int inexact;
    TMP_INIT;

    FLINT_ASSERT(n >= 0 && an >= 1);

    TMP_START;
    N = TMP_ALLOC((nn + dn + (nn - dn + 1) + dn) * sizeof(ulong));
    D = N + nn;
    Q = D + dn;
    R = Q + (nn - dn + 1);

    flint_mpn_zero(N, nn - 1);
    N[nn - 1] = 1;
    D[an] = mpn_lshift(D, a, an, 2 - t);
    FLINT_ASSERT(D[an] != 0);

    flint_mpn_tdiv_qr(Q, R, N, nn, D, dn);
    FLINT_ASSERT(Q[xn] == 0);
    inexact = !flint_mpn_zero_p(R, dn);
    inexact |= (flint_mpn_sqrtrem(y, NULL, Q, xn) != 0);
    TMP_END;
    return inexact;
}
