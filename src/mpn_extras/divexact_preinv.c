/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Exact division by a fixed divisor with a precomputed 2-adic inverse.

    The divisor b = 2^v B^k b' (b' odd, bn' limbs) is stored together with
    binv = b'^(-1) mod B^(bn'). Each division then costs one Hensel block
    division with the classical (block) algorithm: for quotients of at most
    bn' limbs a single low product a' binv mod B^n, and (n / bn') blocks of
    one low and one high bn' x bn' product otherwise; no inverse is ever
    recomputed, which is the main cost of a one-off exact division.
*/

void
flint_mpn_divexact_preinv_init(flint_mpn_divexact_preinv_t pre, mp_srcptr b, mp_size_t bn)
{
    mp_size_t k = 0;

    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    while (b[k] == 0)
        k++;

    pre->k = k;
    pre->bn_orig = bn;
    pre->v = flint_ctz(b[k]);

    b += k;
    bn -= k;

    pre->b = flint_malloc((2 * bn + 1) * sizeof(mp_limb_t));
    pre->binv = pre->b + bn;

    /* the top limb of the odd part may be zero; it is kept so that
       balanced divisions always have n <= bn (a single low product) */
    if (pre->v != 0)
        mpn_rshift(pre->b, b, bn, pre->v);
    else
    {
        flint_mpn_copyi(pre->b, b, bn);
    }

    pre->bn = bn;

    /* one limb beyond bn, so that quotients of bn + 1 limbs (the generic
       balanced case) are still a single low product */
    flint_mpn_binv(pre->binv, pre->b, bn, bn + 1);
}

void
flint_mpn_divexact_preinv_clear(flint_mpn_divexact_preinv_t pre)
{
    flint_free(pre->b);
}

/* top limb of the original divisor: the odd part b' is b / (2^v B^k), so
   b's top limb is recovered from the top two limbs of b' shifted back */
static mp_limb_t
b_top_limb(const flint_mpn_divexact_preinv_t pre)
{
    mp_size_t bn = pre->bn;
    if (pre->v == 0)
        return pre->b[bn - 1];
    return (pre->b[bn - 1] << pre->v) | ((bn > 1) ? (pre->b[bn - 2] >> (FLINT_BITS - pre->v)) : 0);
}

/* q = a / b (an - bn_orig + 1 limbs), assuming b divides a */
void
flint_mpn_divexact_preinv(mp_ptr q, mp_srcptr a, mp_size_t an, const flint_mpn_divexact_preinv_t pre)
{
    mp_size_t n = an - pre->bn_orig + 1, k = pre->k, bn = pre->bn, need;
    mp_ptr as = NULL;
    int truncated = 0;
    TMP_INIT;

    FLINT_ASSERT(an >= pre->bn_orig);
    FLINT_ASSERT(flint_mpn_zero_p(a, k));

    /* the quotient has fewer than n limbs when the top limb of a is below
       that of b (top limbs of the original operands, b + k being the
       original divisor above its zero limbs) */
    if (n > 1 && a[an - 1] < b_top_limb(pre))
    {
        q[n - 1] = 0;
        n--;
    }

    a += k;
    an -= k;

#if FLINT_HAVE_NATIVE_mpn_divexact
    if (bn >= 2 && bn <= 4 && n > bn + 1)
    {
        /* long quotients by tiny divisors: GMP's schoolbook Hensel division
           with a limb inverse beats the block algorithm here; it is given
           the exactly divisible (fully shifted) operands and writes
           an - bn + 1 >= n limbs, the extra ones being zero */
        mp_ptr qq;
        TMP_START;
        as = TMP_ALLOC((an + (an - bn + 1)) * sizeof(mp_limb_t));
        qq = as + an;
        if (pre->v != 0)
            mpn_rshift(as, a, an, pre->v);
        else
            flint_mpn_copyi(as, a, an);
        mpn_divexact(qq, as, an, pre->b, bn);
        flint_mpn_copyi(q, qq, n);
        TMP_END;
        return;
    }
#endif

    /* only the low n limbs of a enter a single low product, and the low
       n + bn limbs a block division */
    need = (n <= bn + 1) ? n : n + bn;
    if (an > need + 1)
    {
        an = need + 1;
        truncated = 1;
    }

    TMP_START;
    if (pre->v != 0)
    {
        as = TMP_ALLOC(an * sizeof(mp_limb_t));
        mpn_rshift(as, a, an, pre->v);
        a = as;
        /* the top shifted limb is incomplete if a was truncated; it is
           never needed */
        if (truncated)
            an--;
    }

    if (bn == 1)
    {
        /* q = a binv mod B^n limb by limb, as in flint_mpn_bdiv_qr_1 */
        mp_limb_t b0 = pre->b[0], binv = pre->binv[0], cy = 0, hi, lo, qi;
        mp_size_t i;

        for (i = 0; i < n; i++)
        {
            qi = (((i < an) ? a[i] : 0) - cy) * binv;
            q[i] = qi;
            umul_ppmm(hi, lo, qi, b0);
            add_ssaaaa(hi, lo, hi, lo, 0, cy);
            cy = hi;
        }
    }
    else if (n <= bn + 1)
    {
        /* a single low product with the (bn + 1)-limb inverse */
        if (an >= n)
        {
            flint_mpn_mullow_n(q, a, pre->binv, n);
        }
        else
        {
            mp_ptr ap = TMP_ALLOC(n * sizeof(mp_limb_t));
            flint_mpn_copyi(ap, a, an);
            flint_mpn_zero(ap + an, n - an);
            flint_mpn_mullow_n(q, ap, pre->binv, n);
        }
    }
    else
    {
        _flint_mpn_bdiv_qr_classical_preinv(q, NULL, a, an, pre->b, bn, pre->binv, n);
    }

    TMP_END;
}
