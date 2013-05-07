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

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

/* Assumes poly1 and poly2 are not length 0. */
void
_nmod_poly_mul_classical(mp_ptr res, mp_srcptr poly1,
                         long len1, mp_srcptr poly2, long len2, nmod_t mod)
{
    long i;
    long log_len = FLINT_BIT_COUNT(len2);
    long bits = FLINT_BITS - (long) mod.norm;

    if (2 * bits + log_len <= FLINT_BITS)
    {
        /* Set res[i] = poly1[i]*poly2[0] */
        mpn_mul_1(res, poly1, len1, poly2[0]);

        if (len2 != 1)
        {
            /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
            mpn_mul_1(res + len1, poly2 + 1, len2 - 1, poly1[len1 - 1]);

            /* out[i+j] += in1[i]*in2[j] */
            for (i = 0; i < len1 - 1; i++)
                mpn_addmul_1(res + i + 1, poly2 + 1, len2 - 1, poly1[i]);
        }

        /* final reduction */
        _nmod_vec_reduce(res, res, len1 + len2 - 1, mod);
    }
    else
    {
        /* Set res[i] = poly1[i]*poly2[0] */
        _nmod_vec_scalar_mul_nmod(res, poly1, len1, poly2[0], mod);
        if (len2 == 1)
            return;

        /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
        _nmod_vec_scalar_mul_nmod(res + len1, poly2 + 1, len2 - 1,
                             poly1[len1 - 1], mod);

        /* out[i+j] += in1[i]*in2[j] */
        for (i = 0; i < len1 - 1; i++)
            _nmod_vec_scalar_addmul_nmod(res + i + 1, poly2 + 1, len2 - 1,
                                    poly1[i], mod);
    }
}

void
nmod_poly_mul_classical(nmod_poly_t res,
                        const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    long len_out;

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
            _nmod_poly_mul_classical(temp->coeffs, poly1->coeffs,
                                     poly1->length, poly2->coeffs,
                                     poly2->length, poly1->mod);
        else
            _nmod_poly_mul_classical(temp->coeffs, poly2->coeffs,
                                     poly2->length, poly1->coeffs,
                                     poly1->length, poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_classical(res->coeffs, poly1->coeffs, poly1->length,
                                     poly2->coeffs, poly2->length, poly1->mod);
        else
            _nmod_poly_mul_classical(res->coeffs, poly2->coeffs, poly2->length,
                                     poly1->coeffs, poly1->length, poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}
