/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

/* assumes res is initially zeroed */
/* assumes that we can write one (zeroed) limb too much */
/* assumes bits >= FLINT_BITS */
static void
_mpn_mod_poly_bit_pack(mp_ptr res, mp_srcptr x, slong len, mp_bitcnt_t bits, mp_size_t nlimbs)
{
    slong i, l, shift;

    for (i = 0; i < len; i++)
    {
        l = (bits * i) / FLINT_BITS;
        shift = (bits * i) % FLINT_BITS;

        if (shift == 0)
            flint_mpn_copyi(res + l, x + i * nlimbs, nlimbs);
        else
            res[l + nlimbs] = mpn_lshift(res + l, x + i * nlimbs, nlimbs, shift);
    }
}

static void
_mpn_mod_poly_bit_unpack(mp_ptr res, mp_srcptr x, slong len, mp_bitcnt_t bits, mp_size_t nlimbs, gr_ctx_t ctx)
{
    slong i, i1, i2, l1, shift, l2;
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS + 3];
    mp_limb_t mask;
    mp_size_t blimbs, tn;

    blimbs = (bits + FLINT_BITS - 1) / FLINT_BITS;

    if (bits % FLINT_BITS == 0)
        mask = ~UWORD(0);
    else
        mask = ((~UWORD(0)) >> (FLINT_BITS - (bits % FLINT_BITS)));

    for (i = 0; i < len; i++)
    {
        /* read bits i1 <= j < i2 */
        i1 = bits * i;
        i2 = bits * (i + 1);

        /* read limbs l1 <= l < l2 */
        l1 = i1 / FLINT_BITS;
        l2 = (i2 + FLINT_BITS - 1) / FLINT_BITS;

        shift = i1 % FLINT_BITS;

        if (shift == 0)
            flint_mpn_copyi(t, x + l1, l2 - l1);
        else
            mpn_rshift(t, x + l1, l2 - l1, shift);

        /* mask off high bits */
        tn = blimbs;
        t[tn - 1] &= mask;
        MPN_NORM(t, tn);
        mpn_mod_set_mpn(res + i * nlimbs, t, tn, ctx);
    }
}

int
_mpn_mod_poly_mullow_KS(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
{
    slong bits, nbits, nlimbs, limbs1, limbs2;
    mp_ptr arr1, arr2, arr;
    int squaring;

    len1 = FLINT_MIN(len1, len);
    len2 = FLINT_MIN(len2, len);

    squaring = (poly1 == poly2 && len1 == len2);

    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    nbits = MPN_MOD_CTX_MODULUS_BITS(ctx);

    bits = 2 * nbits + FLINT_BIT_COUNT(FLINT_MIN(len1, len2));

    limbs1 = (bits * len1 - 1) / FLINT_BITS + 1;
    limbs2 = (bits * len2 - 1) / FLINT_BITS + 1;

    FLINT_ASSERT(limbs1 >= (bits * (len1 - 1) / FLINT_BITS + nlimbs + 1));
    FLINT_ASSERT(limbs2 >= (bits * (len2 - 1) / FLINT_BITS + nlimbs + 1));

    arr1 = flint_calloc(squaring ? limbs1 : limbs1 + limbs2, sizeof(mp_limb_t));
    arr2 = squaring ? arr1 : arr1 + limbs1;

    _mpn_mod_poly_bit_pack(arr1, poly1, len1, bits, nlimbs);
    if (!squaring)
        _mpn_mod_poly_bit_pack(arr2, poly2, len2, bits, nlimbs);

    arr = flint_malloc((limbs1 + limbs2) * sizeof(mp_limb_t));

    if (squaring)
        flint_mpn_sqr(arr, arr1, limbs1);
    else if (limbs1 >= limbs2)
        flint_mpn_mul(arr, arr1, limbs1, arr2, limbs2);
    else
        flint_mpn_mul(arr, arr2, limbs2, arr1, limbs1);

    _mpn_mod_poly_bit_unpack(res, arr, len, bits, nlimbs, ctx);

    flint_free(arr1);
    flint_free(arr);

    return GR_SUCCESS;
}