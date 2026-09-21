/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fixed.h"

#if FLINT_HAVE_NATIVE_mpn_divrem_2
mp_limb_t __gmpn_divrem_2(mp_ptr, mp_size_t, mp_ptr, mp_size_t, mp_srcptr);
#endif

/* floor(B^n / x) for one limb x, written to (q, n + 1) */
static void
_inv_1(mp_ptr q, mp_limb_t x, mp_size_t n)
{
    mp_limb_t one = 1;

    if (x == 1)
    {
        flint_mpn_zero(q, n);
        q[n] = 1;
        return;
    }

#if FLINT_PREINVERT_LIMB_USE_NATIVE
    /* the hardware division chain is latency bound; GMP's pipelined
       mpn_divrem_1 wins for long quotients */
    if (FLINT_MPN_DIVREM_1_USE_HW(n + 1, x))
    {
        mp_limb_t r = 1;
        mp_size_t i;

        q[n] = 0;
        for (i = n - 1; i >= 0; i--)
            udiv_qrnnd(q[i], r, r, 0, x);
        return;
    }
#endif

    /* 1 / x with n fraction limbs */
    mpn_divrem_1(q, n, &one, 1, x);
}

/* floor(B^n / x) for two limbs x (x[1] != 0), written to (q, n) */
static void
_inv_2(mp_ptr q, mp_srcptr x, mp_size_t n)
{
    mp_limb_t d1, d0, r1, r0;
    unsigned int s;

    /* B^n 2^s / (x 2^s) with d = x 2^s normalized; the top two limbs of
       the numerator are (2^s, 0), which is >= d only for x = B */
    s = flint_clz(x[1]);
    d1 = (s == 0) ? x[1] : (x[1] << s) | (x[0] >> (FLINT_BITS - s));
    d0 = x[0] << s;
    r1 = UWORD(1) << s;
    r0 = 0;

    if (r1 > d1 || (r1 == d1 && r0 >= d0))
    {
        /* x = B */
        flint_mpn_zero(q, n - 1);
        q[n - 1] = 1;
        return;
    }

    q[n - 1] = 0;

#if FLINT_HAVE_NATIVE_mpn_divrem_2
    /* GMP's assembly mpn_divrem_2 for long chains */
    if (n + 1 >= FLINT_MPN_DIV_2_GMP_CUTOFF)
    {
        mp_limb_t np[2], dp[2];

        np[0] = r0;
        np[1] = r1;
        dp[0] = d0;
        dp[1] = d1;
        /* (r1, r0) with n - 1 fraction limbs */
        __gmpn_divrem_2(q, n - 1, np, 2, dp);
        return;
    }
#endif
    {
        mp_size_t i;
        mp_limb_t dinv;

#if FLINT_PREINVERT_LIMB_USE_NATIVE
        if (n + 1 < FLINT_MPN_DIV_2_HW_CUTOFF)
        {
            for (i = n - 2; i >= 0; i--)
                FLINT_MPN_UDIV_QR_3BY2_HW(q[i], r1, r0, r1, r0, 0, d1, d0);
            return;
        }
#endif
        dinv = flint_mpn_preinv1(d1, d0);
        for (i = n - 2; i >= 0; i--)
            FLINT_MPN_UDIV_QR_3BY2(q[i], r1, r0, r1, r0, 0, d1, d0, dinv);
    }
}

/*
    floor(B^n / x) without Newton inversion (see flint_mpn_inv): one- and
    two-limb x by a chain of 2/1 or 3/2 divisions of zero limbs, short x
    and short quotients by the register-based division, otherwise the
    quotient of B^n by x without remainder (flint_mpn_tdiv_q, which only
    reads the top 2 (n - xn) + O(1) limbs of x for short quotients).
*/
void
_flint_mpn_inv_basecase(mp_ptr q, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_ptr N;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(n >= xn);
    FLINT_ASSERT(x[xn - 1] != 0);

    if (xn == 1)
    {
        _inv_1(q, x[0], n);
        return;
    }

    if (xn == 2)
    {
        _inv_2(q, x, n);
        return;
    }

    /* short divisors with short quotients: the register-based division,
       with the numerator on the stack */
    if (xn <= FLINT_MPN_DIV_SMALL_BN && n < 3 * FLINT_MPN_DIV_SMALL_BN
        && FLINT_MPN_DIV_SMALL_SHAPE(n + 1, xn, 0))
    {
        mp_limb_t M[3 * FLINT_MPN_DIV_SMALL_BN];

        flint_mpn_zero(M, n);
        M[n] = 1;
        _flint_mpn_tdiv_qr_small(q, NULL, M, n + 1, x, xn);
        return;
    }

    TMP_START;
    N = TMP_ALLOC((n + 1) * sizeof(mp_limb_t));
    flint_mpn_zero(N, n);
    N[n] = 1;
    flint_mpn_tdiv_q(q, N, n + 1, x, xn);
    TMP_END;
}

/*
    Correctly truncated reciprocal: q = floor(B^n / x) for (x, xn) with
    x[xn-1] != 0 and n >= xn, written to (q, n - xn + 2) (the top limb is
    nonzero only when x is a power of B).

    Below the Newton cutoff (in the smaller of xn and the quotient length),
    _flint_mpn_inv_basecase. Newton: viewing x as a fixed-point number
    alpha in [1/B, 1) with xn fraction limbs, B^n / x = (1/alpha) B^(n-xn).
    fixed_inv_newton with p = n - xn + 3 fraction limbs has error at most
    4 B^(-p) / alpha <= 4 B^(-p+1), i.e. 4 B^-2 at the integer scale, so
    the integer part is certified when the first fraction limb lies in
    [2, B-2]; otherwise the quotient is corrected against the explicit
    numerator B^n.
*/
void
flint_mpn_inv(mp_ptr q, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_size_t qn = n - xn + 2, p;
    mp_ptr U, qq;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(n >= xn);
    FLINT_ASSERT(x[xn - 1] != 0);

    /* quotients up to about as long as x use the short division below
       FLINT_MPN_INV_NEWTON_CUTOFF; longer ones have a lower cutoff */
    if (FLINT_MIN(xn, qn) < ((qn <= xn + 2) ? FLINT_MPN_INV_NEWTON_CUTOFF
                                            : FLINT_MPN_INV_NEWTON_LONG_CUTOFF))
    {
        _flint_mpn_inv_basecase(q, x, xn, n);
        return;
    }

    TMP_START;

    p = n - xn + 3;
    U = TMP_ALLOC((p + 3) * sizeof(mp_limb_t));
    U[p + 2] = 0;

    fixed_inv_newton(U, x, xn, p);

    /* integer limb 0 sits at U[p - (n - xn)] = U[3] */
    qq = U + 3;

    if (qq[-1] > 1 && qq[-1] < UWORD_MAX - 1)
    {
        flint_mpn_copyi(q, qq, qn);
    }
    else
    {
        /* verify against B^n: (qq, qn + 1) is within a few units of the
           quotient of the (n + 1)-limb numerator by x */
        mp_ptr N, r;
        N = TMP_ALLOC((2 * (n + 2)) * sizeof(mp_limb_t));
        r = N + n + 2;
        flint_mpn_zero(N, n);
        N[n] = 1;

        if (qq[qn] != 0)
            mpn_sub_1(qq, qq, qn + 1, 1);

        if (qn >= xn)
            flint_mpn_mul(r, qq, qn, x, xn);
        else
            flint_mpn_mul(r, x, xn, qq, qn);

        /* r = qq x has qn + xn = n + 2 limbs */
        while (r[n + 1] != 0 || mpn_cmp(r, N, n + 1) > 0)
        {
            mpn_sub(r, r, n + 2, x, xn);
            mpn_sub_1(qq, qq, qn, 1);
        }

        mpn_sub_n(r, N, r, n + 1);
        r[n + 1] = 0;

        while (!flint_mpn_zero_p(r + xn, n + 2 - xn) || mpn_cmp(r, x, xn) >= 0)
        {
            mpn_add_1(qq, qq, qn, 1);
            mpn_sub(r, r, n + 2, x, xn);
        }

        flint_mpn_copyi(q, qq, qn);
    }

    TMP_END;
}

/*
    Approximate reciprocal with the interface of flint_mpn_inv: sets
    (q, n - xn + 2) to floor(B^n / x) or floor(B^n / x) + 1, i.e.
    flint_mpn_divapprox of B^n by x. The algorithms of flint_mpn_inv without
    its correction steps: the exact division chains for xn <= 2 and the
    register-based division for short x and quotients,
    flint_mpn_divapprox_fraction of B^n (the zero limbs not formed in the
    short and Newton divisions) below the Newton cutoffs of flint_mpn_inv,
    and above, fixed_inv_newton with one guard limb instead of three,
    rounded up by more than its error bound.
*/
void
flint_mpn_invapprox(mp_ptr q, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_size_t qn = n - xn + 2, p;
    mp_ptr U;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(n >= xn);
    FLINT_ASSERT(x[xn - 1] != 0);

    /* the exact division chains (xn <= 2) and register-based division with
       the numerator on the stack are as fast as any approximation */
    if (xn <= 2 || (xn <= FLINT_MPN_DIV_SMALL_BN && n < 3 * FLINT_MPN_DIV_SMALL_BN
        && FLINT_MPN_DIV_SMALL_SHAPE(n + 1, xn, 0)))
    {
        _flint_mpn_inv_basecase(q, x, xn, n);
        return;
    }

    if (FLINT_MIN(xn, qn) < ((qn <= xn + 2) ? FLINT_MPN_INV_NEWTON_CUTOFF
                                            : FLINT_MPN_INV_NEWTON_LONG_CUTOFF))
    {
        mp_limb_t one = 1;
        flint_mpn_divapprox_fraction(q, &one, 1, x, xn, n);
        return;
    }

    /* x as alpha in [1/B, 1) with xn fraction limbs; 1/alpha with
       p = n - xn + 2 fraction limbs has error at most 4 B^(-p + 1), i.e.
       4 units of the guard limb U[1] at the integer scale (the integer
       part is {U + 2, p}). Adding 5 units and truncating gives
       floor(B^n / x) or floor(B^n / x) + 1. */
    TMP_START;
    p = n - xn + 2;
    U = TMP_ALLOC((p + 2) * sizeof(mp_limb_t));
    fixed_inv_newton(U, x, xn, p);
    mpn_add_1(U + 1, U + 1, p + 1, 5);
    flint_mpn_copyi(q, U + 2, qn);
    TMP_END;
}
