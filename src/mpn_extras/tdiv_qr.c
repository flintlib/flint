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

/*
    Euclidean division with quotient and remainder, ported from
    radix_divrem / radix_divrem_newton_karp_markstein / radix_divrem_preinv
    with B = 2^64.

    For (a, an) and (b, bn) with an >= bn >= 1 and b[bn-1] != 0, the
    quotient q = floor(a/b) has n = an - bn + 1 limbs and the remainder
    r = a - q b has bn limbs.

    Viewing a as a fixed-point number in [0, 1) with an fraction limbs and b
    as a fixed-point number in [1/B, 1) with bn fraction limbs, the integer
    quotient is floor((a/b) B^(an-bn)). fixed_div_newton (Karp-Markstein
    division) gives an approximation of a/b with n2 = n + 2 fraction limbs and
    absolute error at most 4 B^(-n2) / b <= 4 B^(-n2+1), i.e. at most 4 B^-2
    at the integer scale. Hence, if the first fraction limb q[-1] of the
    approximation lies in [2, B-2], the integer part is certified and the
    remainder follows from a low product; otherwise the quotient is verified
    and corrected by O(1) steps using a full product.
*/

#define PREINV2_EXTRA_LIMBS 2

/* Cutoffs (tuned on a Skylake-class Xeon with fft_small enabled; see
   profile/p-tdiv_qr.c). Newton division beats GMP for roughly balanced
   operands once both the divisor and the quotient exceed
   FLINT_MPN_TDIV_QR_NEWTON_CUTOFF limbs, where FLINT's multiplication has
   switched to fft_small; block division with a shared approximate inverse
   beats GMP much earlier for unbalanced shapes. */
#ifndef FLINT_MPN_TDIV_QR_UNBALANCED4_CUTOFF   /* an >= 4 bn */
#define FLINT_MPN_TDIV_QR_UNBALANCED4_CUTOFF 64
#endif
#ifndef FLINT_MPN_TDIV_QR_UNBALANCED3_CUTOFF   /* an >= 3 bn */
#define FLINT_MPN_TDIV_QR_UNBALANCED3_CUTOFF 512
#endif

/* Given the approximate quotient q (n + 1 limbs, the top one 0 or 1, and
   the true quotient within a few units), write the exact quotient to Q and
   (if R != NULL) the remainder to R. */
static void
_flint_mpn_tdiv_qr_adjust(mp_ptr Q, mp_ptr R, mp_ptr q, mp_srcptr A,
    mp_size_t An, mp_srcptr B, mp_size_t Bn)
{
    mp_size_t n = An - Bn + 1;
    mp_ptr r;
    TMP_INIT;

    TMP_START;
    r = TMP_ALLOC((An + 1) * sizeof(mp_limb_t));

    if (q[n] != 0)
        mpn_sub_1(q, q, n + 1, 1);

    /* r = q b, An + 1 limbs */
    if (n >= Bn)
        flint_mpn_mul(r, q, n, B, Bn);
    else
        flint_mpn_mul(r, B, Bn, q, n);

    while (r[An] != 0 || mpn_cmp(r, A, An) > 0)
    {
        mpn_sub(r, r, An + 1, B, Bn);
        mpn_sub_1(q, q, n, 1);
    }

    /* r = a - q b (nonnegative) */
    mpn_sub_n(r, A, r, An);
    r[An] = 0;

    while (!flint_mpn_zero_p(r + Bn, An + 1 - Bn) || mpn_cmp(r, B, Bn) >= 0)
    {
        mpn_add_1(q, q, n, 1);
        mpn_sub(r, r, An + 1, B, Bn);
    }

    flint_mpn_copyi(Q, q, n);
    if (R != NULL)
        flint_mpn_copyi(R, r, Bn);

    TMP_END;
}

/* Certified quotient q (n limbs): write it to Q, and if R != NULL the
   remainder a - q b, which is determined by its low Bn limbs. */
static void
_flint_mpn_tdiv_qr_certified(mp_ptr Q, mp_ptr R, mp_srcptr q, mp_srcptr A,
    mp_size_t An, mp_srcptr B, mp_size_t Bn)
{
    mp_size_t n = An - Bn + 1;

    if (R != NULL)
    {
        mp_ptr t;
        TMP_INIT;
        TMP_START;
        t = TMP_ALLOC(Bn * sizeof(mp_limb_t));
        flint_mpn_mulmid(t, q, FLINT_MIN(Bn, n), B, Bn, 0, Bn);
        mpn_sub_n(R, A, t, Bn);
        TMP_END;
    }

    flint_mpn_copyi(Q, q, n);
}

/* the quotient vanishes iff a < b */
static int
_flint_mpn_tdiv_qr_trivial(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    if (flint_mpn_zero_p(A + Bn, An - Bn) && mpn_cmp(A, B, Bn) < 0)
    {
        if (R != NULL)
            flint_mpn_copyi(R, A, Bn);
        flint_mpn_zero(Q, An - Bn + 1);
        return 1;
    }
    return 0;
}

/*
    Division with a precomputed approximate inverse: (Binv, Binvn + 2 limbs)
    is the output of fixed_inv_newton(Binv, B, Bn, Binvn), i.e. an
    approximation of B^Bn / b with Binvn fraction limbs and two integral
    limbs. Requires Binvn >= n + 2 and An >= n + 2 where n = An - Bn + 1.

    With B' the top n2 = n + 2 fraction limbs of Binv (relative error
    e1 <= 5 B^(-n2)), A' the top n2 limbs of a (truncation error
    e2 <= B^(-n2)) and e3 <= n2 B^(-n2+1) the deficit of the middle product,
    the computed approximation differs from a/b by less than
    (n2 + 7) B^(-n2+1), which is far below one unit of the first fraction
    limb q[-1] for B = 2^64.
*/
void
_flint_mpn_tdiv_qr_preinv(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn, mp_srcptr Binv, mp_size_t Binvn)
{
    mp_ptr U, q;
    mp_size_t n, n2;
    TMP_INIT;

    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 1);
    FLINT_ASSERT(B[Bn - 1] != 0);

    if (_flint_mpn_tdiv_qr_trivial(Q, R, A, An, B, Bn))
        return;

    n = An - Bn + 1;
    n2 = n + PREINV2_EXTRA_LIMBS;

    FLINT_ASSERT(Binvn >= n2);
    FLINT_ASSERT(An >= n2);

    TMP_START;
    U = TMP_ALLOC((n2 + 2) * sizeof(mp_limb_t));

    /* n2 fraction limbs x (n2 fraction limbs + 2 integral limbs), keeping
       the top n2 + 2 limbs of the product */
    flint_mpn_mulmid(U, A + An - n2, n2, Binv + Binvn - n2, n2 + 2, n2, 2 * n2 + 2);

    q = U + n2 + 2 - (n + 1);

    if (q[-1] > 1 && q[-1] < UWORD_MAX - 1)
        _flint_mpn_tdiv_qr_certified(Q, R, q, A, An, B, Bn);
    else
        _flint_mpn_tdiv_qr_adjust(Q, R, q, A, An, B, Bn);

    TMP_END;
}

/* Karp-Markstein Newton division */
void
_flint_mpn_tdiv_qr_newton(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    mp_ptr U, q;
    mp_size_t n, n2;
    TMP_INIT;

    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 1);
    FLINT_ASSERT(B[Bn - 1] != 0);

    if (_flint_mpn_tdiv_qr_trivial(Q, R, A, An, B, Bn))
        return;

    n = An - Bn + 1;
    n2 = n + PREINV2_EXTRA_LIMBS;

    TMP_START;
    U = TMP_ALLOC((n2 + 2) * sizeof(mp_limb_t));

    fixed_div_newton(U, A, An, B, Bn, n2);
    FLINT_ASSERT(U[n2 + 1] == 0);

    q = U + n2 + 2 - (n + 1);

    if (q[-1] > 1 && q[-1] < UWORD_MAX - 1)
        _flint_mpn_tdiv_qr_certified(Q, R, q, A, An, B, Bn);
    else
        _flint_mpn_tdiv_qr_adjust(Q, R, q, A, An, B, Bn);

    TMP_END;
}

/*
    Unbalanced division (An > 2 Bn) as a sequence of (2 Bn) x Bn block
    divisions from the top down, sharing one approximate inverse of b with
    Bn + 3 fraction limbs (enough for _flint_mpn_tdiv_qr_preinv on every
    block, whose quotients have at most Bn + 1 limbs). Requires Bn >= 3.
*/
void
_flint_mpn_tdiv_qr_unbalanced(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    mp_ptr Binv, Rt, T;
    mp_size_t Binvn, i, antop;
    mp_limb_t cy;
    TMP_INIT;

    FLINT_ASSERT(An > 2 * Bn);
    FLINT_ASSERT(Bn >= 3);
    FLINT_ASSERT(B[Bn - 1] != 0);

    Binvn = Bn + 3;

    TMP_START;
    Binv = TMP_ALLOC((Binvn + 2 + 3 * Bn) * sizeof(mp_limb_t));
    Rt = Binv + Binvn + 2;
    T = Rt + Bn;

    fixed_inv_newton(Binv, B, Bn, Binvn);

    i = (An + Bn - 1) / Bn - 2;
    antop = An - i * Bn;
    FLINT_ASSERT(antop > Bn && antop <= 2 * Bn);

    _flint_mpn_tdiv_qr_preinv(Q + i * Bn, Rt, A + i * Bn, antop, B, Bn, Binv, Binvn);
    cy = Q[i * Bn];

    for (i--; i >= 0; i--)
    {
        flint_mpn_copyi(T, A + i * Bn, Bn);
        flint_mpn_copyi(T + Bn, Rt, Bn);
        _flint_mpn_tdiv_qr_preinv(Q + i * Bn, Rt, T, 2 * Bn, B, Bn, Binv, Binvn);

        /* the block writes Bn + 1 quotient limbs, the top one being zero;
           restore the limb belonging to the block above */
        Q[(i + 1) * Bn] = cy;
        cy = Q[i * Bn];
    }

    if (R != NULL)
        flint_mpn_copyi(R, Rt, Bn);

    TMP_END;
}

/*
    Division with FLINT's precomputed-inverse routines (flint_mpn_preinvn /
    flint_mpn_divrem_preinvn, which require a normalised divisor): the
    dividend and divisor are shifted left to normalise the divisor and the
    remainder is shifted back. For short divisors (a few dozen limbs) and
    long dividends this beats GMP's mpn_tdiv_qr by up to 40%.
*/
void
_flint_mpn_tdiv_qr_preinvn(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    unsigned int norm;
    mp_ptr as, bs, dinv, rr, qq;
    mp_size_t m = An + 1;
    TMP_INIT;

    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 1);
    FLINT_ASSERT(B[Bn - 1] != 0);

    TMP_START;
    as = TMP_ALLOC((m + 2 * Bn + m + (m + 1 - Bn)) * sizeof(mp_limb_t));
    bs = as + m;
    dinv = bs + Bn;
    rr = dinv + Bn;
    qq = rr + m;

    norm = flint_clz(B[Bn - 1]);
    if (norm != 0)
    {
        as[An] = mpn_lshift(as, A, An, norm);
        mpn_lshift(bs, B, Bn, norm);
    }
    else
    {
        flint_mpn_copyi(as, A, An);
        as[An] = 0;
        flint_mpn_copyi(bs, B, Bn);
    }

    flint_mpn_preinvn(dinv, bs, Bn);
    /* as has m = An + 1 limbs (top limb possibly zero); the quotient has
       m - Bn = An - Bn + 1 limbs and the remainder sits in the low Bn limbs
       of rr; the returned flag is zero since the top limb of as is small */
    flint_mpn_divrem_preinvn(qq, rr, as, m, bs, Bn, dinv);

    flint_mpn_copyi(Q, qq, An - Bn + 1);
    if (R != NULL)
    {
        if (norm != 0)
            mpn_rshift(R, rr, Bn, norm);
        else
            flint_mpn_copyi(R, rr, Bn);
    }

    TMP_END;
}

void
_flint_mpn_tdiv_qr_gmp(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    if (R != NULL)
        mpn_tdiv_qr(Q, R, 0, A, An, B, Bn);
    else
        mpn_tdiv_q(Q, A, An, B, Bn);    /* GMP's mpn_div_q, ~30% cheaper */
}

/* R may be NULL */
void
_flint_mpn_tdiv_qr(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    mp_size_t n = An - Bn + 1;

    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 1);
    FLINT_ASSERT(B[Bn - 1] != 0);

    if (Bn >= FLINT_MPN_TDIV_QR_NEWTON_CUTOFF && n >= FLINT_MPN_TDIV_QR_NEWTON_CUTOFF)
    {
        if (An >= 3 * Bn)
            _flint_mpn_tdiv_qr_unbalanced(Q, R, A, An, B, Bn);
        else
            _flint_mpn_tdiv_qr_newton(Q, R, A, An, B, Bn);
    }
    else if ((Bn >= FLINT_MPN_TDIV_QR_UNBALANCED4_CUTOFF && An >= 4 * Bn)
          || (Bn >= FLINT_MPN_TDIV_QR_UNBALANCED3_CUTOFF && An >= 3 * Bn))
    {
        _flint_mpn_tdiv_qr_unbalanced(Q, R, A, An, B, Bn);
    }
    else if ((Bn >= 32 && An >= 4 * Bn) || (Bn >= 4 && An >= 32 * Bn))
    {
        _flint_mpn_tdiv_qr_preinvn(Q, R, A, An, B, Bn);
    }
    else
    {
        _flint_mpn_tdiv_qr_gmp(Q, R, A, An, B, Bn);
    }
}
