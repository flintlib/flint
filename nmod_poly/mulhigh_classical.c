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
#include "ulong_extras.h"

/* Assumes poly1 and poly2 are not length 0. */
void
_nmod_poly_mulhigh_classical(mp_ptr res, mp_srcptr poly1,
                             len_t len1, mp_srcptr poly2, len_t len2, len_t start,
                             nmod_t mod)
{
    len_t m, n;

    _nmod_vec_zero(res, start);

    if (len1 == 1)              /* Special case if the length of both inputs is 1 */
    {
        if (start == 0)
            res[0] = n_mulmod2_preinv(poly1[0], poly2[0], mod.n, mod.ninv);
    }
    else                        /* Ordinary case */
    {
        len_t i;
        len_t bits = FLINT_BITS - (len_t) mod.norm;
        len_t log_len = FLINT_BIT_COUNT(len2);

        if (2 * bits + log_len <= FLINT_BITS)
        {
            /* Set res[i] = poly1[i]*poly2[0] */
            if (start < len1)
                mpn_mul_1(res + start, poly1 + start, len1 - start, poly2[0]);

            if (len2 != 1)
            {
                /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
                m = FLINT_MAX(len1 - 1, start);

                mpn_mul_1(res + m, poly2 + m - len1 + 1, len2 - 1 + len1 - m,
                          poly1[len1 - 1]);

                /* out[i+j] += in1[i]*in2[j] */
                m = FLINT_MAX(start, len2 - 1);
                for (i = m - len2 + 1; i < len1 - 1; i++)
                {
                    n = FLINT_MAX(i + 1, start);
                    mpn_addmul_1(res + n, poly2 + n - i, len2 + i - n,
                                 poly1[i]);
                }
            }

            _nmod_vec_reduce(res, res, len1 + len2 - 1, mod);
        }
        else
        {
            /* Set res[i] = poly1[i]*poly2[0] */
            if (start < len1)
                _nmod_vec_scalar_mul_nmod(res + start, poly1 + start, len1 - start,
                                     poly2[0], mod);

            if (len2 == 1)
                return;

            /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
            m = FLINT_MAX(len1 - 1, start);
            _nmod_vec_scalar_mul_nmod(res + m, poly2 + m - len1 + 1,
                                 len2 - 1 + len1 - m, poly1[len1 - 1], mod);

            /* out[i+j] += in1[i]*in2[j] */
            m = FLINT_MAX(start, len2 - 1);
            for (i = m - len2 + 1; i < len1 - 1; i++)
            {
                n = FLINT_MAX(i + 1, start);
                _nmod_vec_scalar_addmul_nmod(res + n, poly2 + n - i, len2 + i - n,
                                        poly1[i], mod);
            }
        }
    }
}

void
nmod_poly_mulhigh_classical(nmod_poly_t res,
                            const nmod_poly_t poly1, const nmod_poly_t poly2,
                            len_t start)
{
    len_t len_out = poly1->length + poly2->length - 1;

    if (poly1->length == 0 || poly2->length == 0 || start >= len_out)
    {
        nmod_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mulhigh_classical(temp->coeffs, poly1->coeffs,
                                         poly1->length, poly2->coeffs,
                                         poly2->length, start, poly1->mod);
        else
            _nmod_poly_mulhigh_classical(temp->coeffs, poly2->coeffs,
                                         poly2->length, poly1->coeffs,
                                         poly1->length, start, poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mulhigh_classical(res->coeffs, poly1->coeffs,
                                         poly1->length, poly2->coeffs,
                                         poly2->length, start, poly1->mod);
        else
            _nmod_poly_mulhigh_classical(res->coeffs, poly2->coeffs,
                                         poly2->length, poly1->coeffs,
                                         poly1->length, start, poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}
