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

    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_mul_KS(mp_ptr out, mp_srcptr in1, len_t len1,
                  mp_srcptr in2, len_t len2, mp_bitcnt_t bits, nmod_t mod)
{
    len_t len_out = len1 + len2 - 1, limbs1, limbs2;
    mp_ptr mpn1, mpn2, res;

    if (bits == 0)
    {
        mp_bitcnt_t bits1, bits2, loglen;
        bits1  = _nmod_vec_max_bits(in1, len1);
        bits2  = (in1 == in2) ? bits1 : _nmod_vec_max_bits(in2, len2);
        loglen = FLINT_BIT_COUNT(len2);
        
        bits = bits1 + bits2 + loglen;
    }

    limbs1 = (len1 * bits - 1) / FLINT_BITS + 1;
    limbs2 = (len2 * bits - 1) / FLINT_BITS + 1;

    mpn1 = (mp_ptr) flint_malloc(sizeof(mp_limb_t) * limbs1);
    mpn2 = (in1 == in2) ? mpn1 : (mp_ptr) flint_malloc(sizeof(mp_limb_t) * limbs2);

    _nmod_poly_bit_pack(mpn1, in1, len1, bits);
    if (in1 != in2)
        _nmod_poly_bit_pack(mpn2, in2, len2, bits);

    res = (mp_ptr) flint_malloc(sizeof(mp_limb_t) * (limbs1 + limbs2));

    if (in1 != in2)
        mpn_mul(res, mpn1, limbs1, mpn2, limbs2);
    else
        mpn_mul_n(res, mpn1, mpn1, limbs1);

    _nmod_poly_bit_unpack(out, len_out, res, bits, mod);
    
    flint_free(mpn2);
    if (in1 != in2)
        flint_free(mpn1);

    flint_free(res);
}

void
nmod_poly_mul_KS(nmod_poly_t res,
                 const nmod_poly_t poly1, const nmod_poly_t poly2,
                 mp_bitcnt_t bits)
{
    len_t len_out;

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
