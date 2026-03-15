/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_mulmid_KS(fmpz * res, const fmpz * poly1, slong len1,
                                 const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    int neg1, neg2;
    slong limbs1, limbs2, full_len, loglen;
    slong bits1, bits2, bits;
    ulong *arr1, *arr2, *arr3;
    slong sign = 0;
    slong zero_high = 0;

    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nhi != 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

    /* Low truncation of inputs */
    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    /* High truncation of inputs */
    slong nlo2 = (len1 + len2 - 1) - nlo;

    if (len1 > nlo2)
    {
        slong trunc = len1 - nlo2;
        poly1 += trunc;
        len1 -= trunc;
        nlo -= trunc;
        nhi -= trunc;
    }

    if (len2 > nlo2)
    {
        slong trunc = len2 - nlo2;
        poly2 += trunc;
        len2 -= trunc;
        nlo -= trunc;
        nhi -= trunc;
    }

    FMPZ_VEC_NORM(poly1, len1);
    FMPZ_VEC_NORM(poly2, len2);

    if (len1 == 0 || len2 == 0)
    {
        _fmpz_vec_zero(res, nhi - nlo);
        return;
    }

    neg1 = (fmpz_sgn(poly1 + len1 - 1) > 0) ? 0 : -1;
    neg2 = (fmpz_sgn(poly2 + len2 - 1) > 0) ? 0 : -1;

    full_len = len1 + len2 - 1;

    if (nlo >= FLINT_MIN(nhi, full_len))
    {
       _fmpz_vec_zero(res, nhi - nlo);
        return;
    }

    if (nhi > full_len)
    {
        zero_high = nhi - full_len;
        nhi = full_len;
    }

    bits1 = _fmpz_vec_max_bits(poly1, len1);
    if (bits1 < 0)
    {
        sign = 1;
        bits1 = -bits1;
    }

    if (poly1 == poly2 && len1 == len2)
    {
        bits2 = bits1;
    }
    else
    {
        bits2 = _fmpz_vec_max_bits(poly2, len2);
        if (bits2 < 0)
        {
            sign = 1;
            bits2 = -bits2;
        }
    }

    loglen = FLINT_BIT_COUNT(FLINT_MIN(len1, len2));
    bits = bits1 + bits2 + loglen + sign;

    limbs1 = (bits * len1 - 1) / FLINT_BITS + 1;
    limbs2 = (bits * len2 - 1) / FLINT_BITS + 1;

    if (poly1 == poly2 && len1 == len2)
    {
        arr1 = (nn_ptr) flint_calloc(limbs1, sizeof(ulong));
        arr2 = arr1;
        _fmpz_poly_bit_pack(arr1, poly1, len1, bits, neg1);
    }
    else
    {
        arr1 = (nn_ptr) flint_calloc(limbs1 + limbs2, sizeof(ulong));
        arr2 = arr1 + limbs1;
        _fmpz_poly_bit_pack(arr1, poly1, len1, bits, neg1);
        _fmpz_poly_bit_pack(arr2, poly2, len2, bits, neg2);
    }

    arr3 = (nn_ptr) flint_malloc((limbs1 + limbs2) * sizeof(ulong));

    /* todo: mpn middle products */
    if (arr1 == arr2 && limbs1 == limbs2)
        flint_mpn_sqr(arr3, arr1, limbs1);
    else if (limbs1 > limbs2)
        flint_mpn_mul(arr3, arr1, limbs1, arr2, limbs2);
    else
        flint_mpn_mul(arr3, arr2, limbs2, arr1, limbs1);

    if (sign)
        _fmpz_poly_bit_unpack(res, nlo, nhi, arr3, bits, neg1 ^ neg2);
    else
        _fmpz_poly_bit_unpack_unsigned(res, nlo, nhi, arr3, bits);

    if (zero_high != 0)
        _fmpz_vec_zero(res + nhi - nlo, zero_high);

    flint_free(arr1);
    flint_free(arr3);
}

void
fmpz_poly_mulmid_KS(fmpz_poly_t res, const fmpz_poly_t poly1,
                                            const fmpz_poly_t poly2, slong nlo, slong nhi)
{
    slong len1 = fmpz_poly_length(poly1);
    slong len2 = fmpz_poly_length(poly2);
    slong len;

    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi >= 0);

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
    {
        fmpz_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;
    fmpz_poly_fit_length(res, len);
    _fmpz_poly_mulmid_KS(res->coeffs, poly1->coeffs, poly1->length,
                                poly2->coeffs, poly2->length, nlo, nhi);
    _fmpz_poly_set_length(res, len);
    _fmpz_poly_normalise(res);
}

