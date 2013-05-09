/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_mul_KS(fmpz * res, const fmpz * poly1, len_t len1,
                              const fmpz * poly2, len_t len2)
{
    const len_t in1_len = len1, in2_len = len2;
    int neg1, neg2;
    len_t limbs1, limbs2, loglen;
    len_t bits1, bits2, bits;
    mp_limb_t *arr1, *arr2, *arr3;
    len_t sign = 0;

    FMPZ_VEC_NORM(poly1, len1);
    FMPZ_VEC_NORM(poly2, len2);

    if (!len1 | !len2)
    {
        if (in1_len + in2_len - 1 > 0)
            _fmpz_vec_zero(res, in1_len + in2_len - 1);
        return;
    }

    neg1 = (fmpz_sgn(poly1 + len1 - 1) > 0) ? 0 : -1;
    neg2 = (fmpz_sgn(poly2 + len2 - 1) > 0) ? 0 : -1;

    bits1 = _fmpz_vec_max_bits(poly1, len1);
    if (bits1 < 0)
    {
        sign = 1;
        bits1 = -bits1;
    }

    if (poly1 != poly2)
    {
        bits2 = _fmpz_vec_max_bits(poly2, len2);
        if (bits2 < 0)
        {
            sign = 1;
            bits2 = -bits2;
        }
    }
    else
        bits2 = bits1;

    loglen = FLINT_BIT_COUNT(FLINT_MIN(len1, len2));
    bits = bits1 + bits2 + loglen + sign;

    limbs1 = (bits * len1 - 1) / FLINT_BITS + 1;
    limbs2 = (bits * len2 - 1) / FLINT_BITS + 1;

    if (poly1 == poly2)
    {
        arr1 = (mp_limb_t *) flint_calloc(limbs1, sizeof(mp_limb_t));
        arr2 = arr1;
        _fmpz_poly_bit_pack(arr1, poly1, len1, bits, neg1);
    }
    else
    {
        arr1 = (mp_limb_t *) flint_calloc(limbs1 + limbs2, sizeof(mp_limb_t));
        arr2 = arr1 + limbs1;
        _fmpz_poly_bit_pack(arr1, poly1, len1, bits, neg1);
        _fmpz_poly_bit_pack(arr2, poly2, len2, bits, neg2);
    }

    arr3 = (mp_limb_t *) flint_malloc((limbs1 + limbs2) * sizeof(mp_limb_t));

    if (limbs1 == limbs2)
        mpn_mul_n(arr3, arr1, arr2, limbs1);
    else if (limbs1 > limbs2)
        mpn_mul(arr3, arr1, limbs1, arr2, limbs2);
    else
        mpn_mul(arr3, arr2, limbs2, arr1, limbs1);

    if (sign)
        _fmpz_poly_bit_unpack(res, len1 + len2 - 1, arr3, bits, neg1 ^ neg2);
    else
        _fmpz_poly_bit_unpack_unsigned(res, len1 + len2 - 1, arr3, bits);

    if ((len1 < in1_len) | (len2 < in2_len))
        _fmpz_vec_zero(res + (len1 + len2 - 1), (in1_len - len1) + (in2_len - len2));

    flint_free(arr1);
    flint_free(arr3);
}

void
fmpz_poly_mul_KS(fmpz_poly_t res,
                 const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const len_t len1 = poly1->length;
    const len_t len2 = poly2->length;
    const len_t rlen = len1 + len2 - 1;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        fmpz_poly_fit_length(res, rlen);
        if (len1 >= len2)
            _fmpz_poly_mul_KS(res->coeffs, poly1->coeffs, len1,
                              poly2->coeffs, len2);
        else
            _fmpz_poly_mul_KS(res->coeffs, poly2->coeffs, len2,
                              poly1->coeffs, len1);
        _fmpz_poly_set_length(res, rlen);
    }
}
