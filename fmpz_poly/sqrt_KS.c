/*
    Copyright (C) 2008, 2009, 2018 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#define FMPZ_POLY_SQRT_KS_HEURISTIC_BITS 3

int
_fmpz_poly_sqrt_KS(fmpz *rop, const fmpz *op, slong len)
{
    slong i, len2, m, rlimbs;
    int result = 1;
    slong bits, bits2, limbs, limbs2, loglen;
    mp_limb_t *arr, *arr2, *arr3;

    /* the degree must be even */
    if (len % 2 == 0)
        return 0;

    /* valuation must be even, and then can be reduced to 0 */
    while (fmpz_is_zero(op))
    {
        if (!fmpz_is_zero(op + 1))
            return 0;

        fmpz_zero(rop);
        op += 2;
        len -= 2;
        rop++;
    }

    /* check whether a square root exists modulo 2 */
    m = (len + 1) / 2;
    
    for (i = ((m - 1) | 1); i < len; i += 2)
        if (!fmpz_is_even(op + i))
            return 0;

    for (i = 1; i < ((m - 1) | 1); i += 2)
        if (!fmpz_is_even(op + i))
            return 0;

    /* check endpoints */
    if (!fmpz_is_square(op))
        return 0;

    if (len > 1 && !fmpz_is_square(op + len - 1))
        return 0;

    bits = FLINT_ABS(_fmpz_vec_max_bits(op, len));

    len2 = (len + 1) / 2;

    bits += FMPZ_POLY_SQRT_KS_HEURISTIC_BITS + FLINT_BIT_COUNT(len);

    limbs  = (bits * len - 1) / FLINT_BITS + 1;

    arr = (mp_limb_t *) flint_calloc(limbs, sizeof(mp_limb_t));

    _fmpz_poly_bit_pack(arr, op, len, bits, 0);

    limbs2  = (bits * len2 - 1) / FLINT_BITS + 1;
    arr2 = (mp_limb_t *) flint_calloc(limbs2, sizeof(mp_limb_t));

    arr3 = (mp_limb_t *) flint_calloc(limbs, sizeof(mp_limb_t));

    while (limbs != 0 && arr[limbs - 1] == 0)
        limbs--;

    rlimbs = mpn_sqrtrem(arr2, arr3, arr, limbs);

    loglen = FLINT_BIT_COUNT(len2);
    
    if (rlimbs != 0)
        result = 0;
    else
    {
        _fmpz_poly_bit_unpack(rop, len2, arr2, bits, 0);

        bits2 = _fmpz_vec_max_bits(rop, len2);

        if (2*FLINT_ABS(bits2) + loglen + 1 > bits)
            result = -1;
    }

    flint_free(arr);
    flint_free(arr2);
    flint_free(arr3);

    return result;
}

int
fmpz_poly_sqrt_KS(fmpz_poly_t b, const fmpz_poly_t a)
{
    slong blen, len = a->length;
    int result;

    if (len % 2 == 0)
    {
        fmpz_poly_zero(b);
        return len == 0;
    }

    if (b == a)
    {
        fmpz_poly_t tmp;
        fmpz_poly_init(tmp);
        result = fmpz_poly_sqrt_KS(tmp, a);
        fmpz_poly_swap(b, tmp);
        fmpz_poly_clear(tmp);
        return result;
    }

    blen = len / 2 + 1;
    fmpz_poly_fit_length(b, blen);
    _fmpz_poly_set_length(b, blen);
    result = _fmpz_poly_sqrt_KS(b->coeffs, a->coeffs, len);
    if (result <= 0)
        _fmpz_poly_set_length(b, 0);

    return result;
}

