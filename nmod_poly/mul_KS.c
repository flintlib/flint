/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_mul_KS(mp_ptr out, mp_srcptr in1, slong len1,
                  mp_srcptr in2, slong len2, flint_bitcnt_t bits, nmod_t mod)
{
    slong len_out = len1 + len2 - 1, limbs1, limbs2;
    mp_ptr tmp, mpn1, mpn2, res;
    int squaring;
    TMP_INIT;

    squaring = (in1 == in2 && len1 == len2);

    if (bits == 0)
    {
        flint_bitcnt_t bits1, bits2, loglen;

        /* Look at the actual bits of the input? This slows down the generic
        case. Are there situations where we care enough about special input? */
#if 0
        bits1  = _nmod_vec_max_bits2(in1, len1);
        bits2  = squaring ? bits1 : _nmod_vec_max_bits2(in2, len2);
#else
        bits1 = FLINT_BITS - (slong) mod.norm;
        bits2 = bits1;
#endif

        loglen = FLINT_BIT_COUNT(len2);
        bits = bits1 + bits2 + loglen;
    }

    limbs1 = (len1 * bits - 1) / FLINT_BITS + 1;
    limbs2 = (len2 * bits - 1) / FLINT_BITS + 1;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(mp_limb_t) * (limbs1 + limbs2 + limbs1 + (squaring ? 0 : limbs2)));
    res = tmp;
    mpn1 = tmp + limbs1 + limbs2;
    mpn2 = squaring ? mpn1 : (mpn1 + limbs1);

    _nmod_poly_bit_pack(mpn1, in1, len1, bits);
    if (!squaring)
        _nmod_poly_bit_pack(mpn2, in2, len2, bits);

    if (squaring)
        mpn_sqr(res, mpn1, limbs1);
    else
        mpn_mul(res, mpn1, limbs1, mpn2, limbs2);

    _nmod_poly_bit_unpack(out, len_out, res, bits, mod);

    TMP_END;
}

void
nmod_poly_mul_KS(nmod_poly_t res,
                 const nmod_poly_t poly1, const nmod_poly_t poly2,
                 flint_bitcnt_t bits)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS(temp->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, bits,
                              poly1->mod);
        else
            _nmod_poly_mul_KS(temp->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length, bits,
                              poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS(res->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, bits,
                              poly1->mod);
        else
            _nmod_poly_mul_KS(res->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length, bits,
                              poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}
