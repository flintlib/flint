/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_poly.h"

void
_nmod_poly_mulmid_KS(nn_ptr out, nn_srcptr in1, slong len1,
            nn_srcptr in2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    slong limbs1, limbs2, dlen;
    nn_ptr tmp, mpn1, mpn2, res;
    ulong bits;
    int squaring;
    TMP_INIT;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    if (nlo != 0)
    {
        slong nlo2 = (len1 + len2 - 1) - nlo;

        if (len1 > nlo2)
        {
            slong trunc = len1 - nlo2;
            in1 += trunc;
            len1 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }

        if (len2 > nlo2)
        {
            slong trunc = len2 - nlo2;
            in2 += trunc;
            len2 -= trunc;
            nlo -= trunc;
            nhi -= trunc;
        }
    }

    squaring = (in1 == in2 && len1 == len2);

    dlen = FLINT_MIN(len1, len2);
    bits = 2 * NMOD_BITS(mod) + FLINT_BIT_COUNT(dlen);
    if (bits <= FLINT_BITS)
    {
        ulong m = (mod.n - 1) * (mod.n - 1) * dlen;
        bits = FLINT_BIT_COUNT(m);
        bits = FLINT_MAX(bits, 1);
    }

    limbs1 = (len1 * bits - 1) / FLINT_BITS + 1;
    limbs2 = (len2 * bits - 1) / FLINT_BITS + 1;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(ulong) * (limbs1 + limbs2 + limbs1 + (squaring ? 0 : limbs2)));
    res = tmp;
    mpn1 = tmp + limbs1 + limbs2;
    mpn2 = squaring ? mpn1 : (mpn1 + limbs1);

    _nmod_poly_bit_pack(mpn1, in1, len1, bits);
    if (!squaring)
        _nmod_poly_bit_pack(mpn2, in2, len2, bits);

    if (squaring)
        flint_mpn_sqr(res, mpn1, limbs1);
    else if (limbs1 >= limbs2)
        flint_mpn_mul(res, mpn1, limbs1, mpn2, limbs2);
    else
        flint_mpn_mul(res, mpn2, limbs2, mpn1, limbs1);

    _nmod_poly_bit_unpack(out, nlo, nhi, res, bits, mod);

    TMP_END;
}

void
_nmod_poly_mullow_KS(nn_ptr out, nn_srcptr in1, slong len1,
                  nn_srcptr in2, slong len2, slong n, nmod_t mod)
{
    _nmod_poly_mulmid_KS(out, in1, len1, in2, len2, 0, n, mod);
}

void
_nmod_poly_mul_KS(nn_ptr out, nn_srcptr in1, slong len1,
                  nn_srcptr in2, slong len2, nmod_t mod)
{
    _nmod_poly_mulmid_KS(out, in1, len1, in2, len2, 0, len1 + len2 - 1, mod);
}

void
nmod_poly_mulmid_KS(nmod_poly_t res,
    const nmod_poly_t poly1, const nmod_poly_t poly2, slong nlo, slong nhi)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len;

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
    {
        nmod_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;

    nmod_poly_fit_length(res, len);
    _nmod_poly_mulmid_KS(res->coeffs, poly1->coeffs, len1,
        poly2->coeffs, len2, nlo, nhi, poly1->mod);
    res->length = len;
    _nmod_poly_normalise(res);
}

void
nmod_poly_mullow_KS(nmod_poly_t res,
    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    len = FLINT_MIN(n, len1 + len2 - 1);
    nmod_poly_fit_length(res, len);
    _nmod_poly_mullow_KS(res->coeffs, poly1->coeffs, len1,
        poly2->coeffs, len2, len, poly1->mod);
    res->length = len;
    _nmod_poly_normalise(res);
}

void
nmod_poly_mul_KS(nmod_poly_t res,
    const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len;

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    len = len1 + len2 - 1;

    nmod_poly_fit_length(res, len);
    _nmod_poly_mul_KS(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, poly1->mod);
    res->length = len;
    _nmod_poly_normalise(res);
}

