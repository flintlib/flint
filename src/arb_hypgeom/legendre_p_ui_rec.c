/*
    Copyright (C) 2017, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "arb_hypgeom.h"

/*
    The three-term recurrence scaled by its denominators: with
    u_k = k! P_k(x), (k+1) P_(k+1) = (2k+1) x P_k - k P_(k-1) becomes

        u_(k+1) = (2k+1) x u_k - k^2 u_(k-1),

    with no division.  For 0 <= x < 1 given as an L-limb fraction
    X / beta^L (beta = 2^FLINT_BITS), u_(k-1) and u_k are kept as L-limb
    magnitudes with signs at a common limb exponent S (u = +-M beta^S,
    M < beta^L).  A step is one high product x u_k ~ mulhigh(X, M) beta^S,
    a multiplication by a word and an addmul or submul by a word; when
    the result needs L + 1 limbs, both values drop their lowest limb.  At the end,
    P_n = u_n / n! and P_(n-1) = u_(n-1) / (n-1)!.

    Error bound.  Let p_k = u_k / k! for the computed u_k (after any
    truncation).  Then (k+1) p_(k+1) = (2k+1) x p_k - k p_(k-1) +
    eta_(k+1), where k! eta_(k+1) collects the errors of step k: the high
    product (less than 2 units beta^S_k [AJ2025], times 2k+1) and, if the
    step drops a limb, the truncations of u_(k+1) and of u_k (less than
    one unit beta^S' of the new exponent each, the latter times
    (2k+1) x <= 2k+1, since u_k was used before its truncation).

    The exponents satisfy beta^S_k <= beta^(1-L) k! (1 + d) during step
    k, and beta^S' <= beta^(1-L) (k+1)! (1 + d) after a drop at step k,
    since a drop happens only when |u_(k+1)| >= beta^(L - 1 + S') and
    |p_j| <= 1 + d, with d the largest error |p_j - P_j(x)|, j <= k + 1
    (initially beta^S = beta^(1-L)).  Hence
        |eta_(k+1)| <= beta^(1-L) (1 + d) (2 (k+1)^2 + 4 (2k+1)),
    and |eta_1| < beta^(1-L) from the initial truncation of u_1.

    By Proposition 5 of [JM2018] (the generating function of the errors
    is p(z) int_0^z eta(w) p(w) dw, with p(z) = sum P_k(x) z^k and
    |P_k(x)| <= 1), |p_n - P_n(x)| <= sum_j |eta_j| (H_n - H_(j-1)), and
    H_n - H_(j-1) <= (n - j + 1) / j.  Summing, with 1 + d <= 2,

        |p_n - P_n(x)|, |p_(n-1) - P_(n-1)(x)|
            <= beta^(1-L) (2 n (n+1) (n+2) / 3 + 8 n^2 + n) = E,

    and d <= 1 holds by induction when E <= 1, which the working
    precision ensures.
*/

/* Runs the recurrence to u_n, u_(n-1) (magnitudes Mn, Mn1 of L limbs,
   negative flags sn, sn1, common limb exponent *e); X is the L-limb
   fraction, L >= 2.  tmp has space for 3 (L + 1) limbs. */
static void
_legendre_rec_mpn(nn_ptr Mn, int * sn, nn_ptr Mn1, int * sn1, slong * e,
    nn_srcptr X, slong L, ulong n, nn_ptr tmp)
{
    nn_ptr T, A, Bv, t;
    ulong k, hi, bw;
    slong S;
    int sa, sb, sw;

    T = tmp;
    A = T + (L + 1);
    Bv = A + (L + 1);

    /* u_0 = 1 = beta^(L-1) beta^S, u_1 = x truncated to a multiple of
       beta^S, S = 1 - L */
    S = 1 - L;
    flint_mpn_zero(A, L);
    A[L - 1] = 1;
    flint_mpn_copyi(Bv, X + 1, L - 1);
    Bv[L - 1] = 0;
    sa = 0;
    sb = 0;

    for (k = 1; k < n; k++)
    {
        /* T = (2k+1) mulhigh(X, |u_k|) */
        flint_mpn_mulhigh_n(T, X, Bv, L);
        T[L] = mpn_mul_1(T, T, L, 2 * k + 1);

        /* u_(k+1) = (2k+1) x u_k - k^2 u_(k-1): the first term has the
           sign of u_k, the second the opposite of that of u_(k-1) */
        sw = sb;
        if (sb != sa)
        {
            T[L] += mpn_addmul_1(T, A, L, k * k);
        }
        else
        {
            hi = T[L];
            bw = mpn_submul_1(T, A, L, k * k);
            T[L] = hi - bw;
            if (hi < bw)
            {
                mpn_neg(T, T, L + 1);
                sw = !sb;
            }
        }

        /* keep L limbs: drop the lowest limb of u_(k+1) and u_k */
        if (T[L] != 0)
        {
            flint_mpn_copyi(T, T + 1, L);
            flint_mpn_copyi(Bv, Bv + 1, L - 1);
            Bv[L - 1] = 0;
            S++;
        }

        /* (u_(k-1), u_k) <- (u_k, u_(k+1)) */
        t = A;
        A = Bv;
        Bv = T;
        T = t;
        sa = sb;
        sb = sw;
    }

    flint_mpn_copyi(Mn, Bv, L);
    flint_mpn_copyi(Mn1, A, L);
    *sn = sb;
    *sn1 = sa;
    *e = S;
}

/* the largest k with k! < 2^FLINT_BITS */
#if FLINT_BITS == 64
#define FAC_WORD_MAX 20
#else
#define FAC_WORD_MAX 12
#endif

/* k! for k <= FAC_WORD_MAX */
static ulong
_fac_word(ulong k)
{
    ulong f = 1, i;
    for (i = 2; i <= k; i++)
        f *= i;
    return f;
}

/* res = (-1)^neg M beta^e / k!, plus the error bound E; f is k! if
   k > FAC_WORD_MAX */
static void
_legendre_rec_finish(arb_t res, nn_srcptr M, slong L, int neg, slong e,
    ulong k, const arb_t f, const mag_t E, slong prec)
{
    fmpz_t z;

    fmpz_init(z);
    fmpz_set_ui_array(z, M, L);
    if (neg)
        fmpz_neg(z, z);
    arb_set_fmpz(res, z);
    arb_mul_2exp_si(res, res, FLINT_BITS * e);
    if (k <= FAC_WORD_MAX)
        arb_div_ui(res, res, _fac_word(k), prec);
    else
        arb_div(res, res, f, prec);
    arb_add_error_mag(res, E);
    fmpz_clear(z);
}

void
arb_hypgeom_legendre_p_ui_rec(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong prec)
{
    slong wp, L, e;
    nn_ptr X, Mn, Mn1, tmp;
    fmpz_t fx;
    int inexact, sn, sn1, xneg;
    arb_t t, u, f, x2sub1;
    mag_t err1, err2, xrad, E;
    TMP_INIT;

    if (n > (UWORD(1) << (FLINT_BITS / 2 - 1)))
    {
        if (res != NULL) arb_indeterminate(res);
        if (res_prime != NULL) arb_indeterminate(res_prime);
        return;
    }

    if (n == 0)
    {
        if (res != NULL) arb_one(res);
        if (res_prime != NULL) arb_zero(res_prime);
        return;
    }

    /* the error analysis assumes |x| < 1 on the whole ball (decided
       exactly: a magnitude bound would round 1 - 2^-k up to 1) */
    {
        arf_t ub;
        int inside;
        arf_init(ub);
        arb_get_abs_ubound_arf(ub, x, ARF_PREC_EXACT);
        inside = arf_cmpabs_2exp_si(ub, 0) < 0;
        arf_clear(ub);
        if (!inside)
        {
            arb_hypgeom_legendre_p_ui_one(res, res_prime, n, x, n + 1, prec);
            return;
        }
    }

    mag_init(xrad);

    fmpz_init(fx);
    arb_init(t);
    arb_init(u);
    arb_init(f);
    arb_init(x2sub1);
    mag_init(err1);
    mag_init(err2);
    mag_init(E);

    /* absolute precision: prec bits relative to |x| for small x (P_n(x)
       is about x times a constant for odd n), and the bits lost to the
       error bound E, about n^3 */
    wp = -arf_abs_bound_lt_2exp_si(arb_midref(x));
    wp = FLINT_MAX(wp, 0);
    wp = FLINT_MIN(wp, prec);
    wp += prec + 3 * FLINT_BIT_COUNT(n) + 4;
    L = (wp + FLINT_BITS - 1) / FLINT_BITS + 1;

    arb_mul(x2sub1, x, x, ARF_PREC_EXACT);
    arb_neg(x2sub1, x2sub1);
    arb_add_ui(x2sub1, x2sub1, 1, wp);

    /* X = trunc(|mid(x)| beta^L); P_k(-x) = (-1)^k P_k(x) */
    xneg = arf_sgn(arb_midref(x)) < 0;
    inexact = arf_get_fmpz_fixed_si(fx, arb_midref(x), -FLINT_BITS * L);
    fmpz_abs(fx, fx);

    mag_set(xrad, arb_radref(x));
    if (inexact)
        mag_add_ui_2exp_si(xrad, xrad, 1, -FLINT_BITS * L);

    TMP_START;
    X = TMP_ALLOC(sizeof(ulong) * (6 * (L + 1)));
    Mn = X + (L + 1);
    Mn1 = Mn + (L + 1);
    tmp = Mn1 + (L + 1);
    fmpz_get_ui_array(X, L, fx);

    _legendre_rec_mpn(Mn, &sn, Mn1, &sn1, &e, X, L, n, tmp);

    /* E = beta^(1-L) (2 n (n+1) (n+2) / 3 + 8 n^2 + n) */
    mag_set_ui(E, n);
    mag_mul_ui(E, E, n + 1);
    mag_mul_ui(E, E, n + 2);
    mag_mul_2exp_si(E, E, 1);
    mag_div_ui(E, E, 3);
    mag_set_ui(err1, n);
    mag_mul_ui(err1, err1, 8 * n + 1);
    mag_add(E, E, err1);
    mag_mul_2exp_si(E, E, FLINT_BITS * (1 - L));

    if (!mag_is_zero(xrad))
    {
        arb_hypgeom_legendre_p_ui_deriv_bound(err1, err2, n, x, x2sub1);
        mag_mul(err1, err1, xrad);
        mag_mul(err2, err2, xrad);
    }
    else
    {
        mag_zero(err1);
        mag_zero(err2);
    }

    /* P_n = u_n / n!, P_(n-1) = u_(n-1) / (n-1)! */
    if (n > FAC_WORD_MAX)
        arb_fac_ui(f, n, wp);
    _legendre_rec_finish(t, Mn, L, sn ^ (xneg && (n % 2 == 1)), e, n, f, E, wp);
    mag_add(arb_radref(t), arb_radref(t), err1);

    if (res_prime != NULL)
    {
        /* P' = n (P[n-1] - x P) / (1 - x^2) */
        if (n - 1 > FAC_WORD_MAX)
            arb_div_ui(f, f, n, wp);
        _legendre_rec_finish(u, Mn1, L, sn1 ^ (xneg && (n % 2 == 0)), e, n - 1, f, E, wp);

        arb_submul(u, t, x, wp);
        arb_div(res_prime, u, x2sub1, wp);
        arb_mul_ui(res_prime, res_prime, n, prec);

        mag_add(arb_radref(res_prime), arb_radref(res_prime), err2);
    }

    if (res != NULL)  /* done last since x may be aliased with res */
    {
        arb_set_round(res, t, prec);
    }

    TMP_END;
    fmpz_clear(fx);
    arb_clear(t);
    arb_clear(u);
    arb_clear(f);
    arb_clear(x2sub1);
    mag_clear(err1);
    mag_clear(err2);
    mag_clear(E);
    mag_clear(xrad);
}
