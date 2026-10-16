/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"

void
_fmpz_poly_mulmid_mpn_bits(fmpz * res, const fmpz * poly1, slong len1, slong bits1,
    const fmpz * poly2, slong len2, slong bits2, slong nlo, slong nhi)
{
    slong len, n1, n2, s, ub1, ub2, alloc;
    int sgn1, sgn2, method, squaring, direct1, direct2;
    nn_ptr a, b, t;
    mpn_poly_mul_params_t P;
    TMP_INIT;

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
    {
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
    }

    len = nhi - nlo;

    if (len1 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly2 + nlo, len, poly1);
        return;
    }

    if (len2 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly1 + nlo, len, poly2);
        return;
    }

    if (bits1 == 0 || bits2 == 0)
    {
        _fmpz_vec_zero(res, len);
        return;
    }

    squaring = (poly1 == poly2 && len1 == len2);

    sgn1 = bits1 < 0;
    sgn2 = bits2 < 0;
    ub1 = FLINT_ABS(bits1);
    ub2 = FLINT_ABS(bits2);

    /* representation, limb counts and output size */
    _flint_mpn_poly_mulmid_params(P, len1, bits1, len2, bits2, nlo, nhi, squaring);
    method = P->method;
    n1 = P->nlimbs1;
    n2 = P->nlimbs2;
    s = P->slimbs;

    /* When all coefficients are small fmpz values, the fmpz arrays are
       already valid arrays of 1-limb two's complement integers. */
    direct1 = (n1 == 1 && ub1 <= SMALL_FMPZ_BITCOUNT_MAX && (method == FLINT_MPN_POLY_MUL_SIGNED || method == FLINT_MPN_POLY_MUL_UNSIGNED));
    direct2 = (n2 == 1 && ub2 <= SMALL_FMPZ_BITCOUNT_MAX && (method == FLINT_MPN_POLY_MUL_SIGNED || method == FLINT_MPN_POLY_MUL_UNSIGNED));

    /* sign-magnitude coefficients have an extra limb */
    alloc = (direct1 ? 0 : len1 * (n1 + (method == FLINT_MPN_POLY_MUL_SIGNMAG)))
            + ((squaring || direct2) ? 0 : len2 * (n2 + (method == FLINT_MPN_POLY_MUL_SIGNMAG))) + len * s;

    TMP_START;
    t = TMP_ALLOC(sizeof(ulong) * alloc);

    if (direct1)
    {
        a = (nn_ptr) poly1;
    }
    else if (method == FLINT_MPN_POLY_MUL_SIGNMAG)
    {
        a = t;
        t += len1 * (n1 + 1);
        _fmpz_vec_get_limbs_signmag(a, poly1, len1, n1);
    }
    else
    {
        a = t;
        t += len1 * n1;
        _fmpz_vec_get_limbs(a, poly1, len1, n1, (method == FLINT_MPN_POLY_MUL_BIAS && sgn1) ? ub1 : -1);
    }

    if (squaring)
    {
        b = a;
    }
    else if (direct2)
    {
        b = (nn_ptr) poly2;
    }
    else if (method == FLINT_MPN_POLY_MUL_SIGNMAG)
    {
        b = t;
        t += len2 * (n2 + 1);
        _fmpz_vec_get_limbs_signmag(b, poly2, len2, n2);
    }
    else
    {
        b = t;
        t += len2 * n2;
        _fmpz_vec_get_limbs(b, poly2, len2, n2, (method == FLINT_MPN_POLY_MUL_BIAS && sgn2) ? ub2 : -1);
    }

    _flint_mpn_poly_mulmid(t, a, len1, n1, P->norm1, b, len2, n2, P->norm2, nlo, nhi, s, method);

    _fmpz_poly_mpn_set_outputs(res, t, s, nlo, nhi, method,
        a, len1, n1, b, len2, n2, ub1, ub2, sgn1, sgn2, squaring);

    TMP_END;
}

void
_fmpz_poly_mulmid_mpn(fmpz * res, const fmpz * poly1, slong len1,
                      const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    slong bits1, bits2;

    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nhi != 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

    /* Truncate inputs before computing the bit bounds */
    {
        slong nlo2;

        len1 = FLINT_MIN(len1, nhi);
        len2 = FLINT_MIN(len2, nhi);
        nlo2 = (len1 + len2 - 1) - nlo;

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
    }

    bits1 = _fmpz_vec_max_bits(poly1, len1);
    bits2 = (poly1 == poly2 && len1 == len2) ? bits1 : _fmpz_vec_max_bits(poly2, len2);

    _fmpz_poly_mulmid_mpn_bits(res, poly1, len1, bits1, poly2, len2, bits2, nlo, nhi);
}

void
_fmpz_poly_mul_mpn(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
{
    _fmpz_poly_mulmid_mpn(res, poly1, len1, poly2, len2, 0, len1 + len2 - 1);
}

void
fmpz_poly_mulmid_mpn(fmpz_poly_t res, const fmpz_poly_t poly1,
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

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, len);
        _fmpz_poly_mulmid_mpn(t->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, nlo, nhi);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, len);
        _fmpz_poly_mulmid_mpn(res->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, nlo, nhi);
    }

    _fmpz_poly_set_length(res, len);
    _fmpz_poly_normalise(res);
}

void
fmpz_poly_mul_mpn(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    fmpz_poly_mulmid_mpn(res, poly1, poly2, 0, fmpz_poly_length(poly1) + fmpz_poly_length(poly2) - 1);
}
