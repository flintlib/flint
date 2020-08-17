/*
    Copyright (C) 2008, 2009 William Hart

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
#include "ulong_extras.h"

/* Assumes poly1 and poly2 are not length 0 and 0 < trunc <= len1 + len2 - 1 */
void
_nmod_poly_mullow_classical(mp_ptr res, mp_srcptr poly1, slong len1,
                            mp_srcptr poly2, slong len2, slong trunc, nmod_t mod)
{
    if (len1 == 1 || trunc == 1)    /* Special case if the length of output is 1 */
    {
        res[0] = n_mulmod2_preinv(poly1[0], poly2[0], mod.n, mod.ninv);
    }
    else                        /* Ordinary case */
    {
        slong i;

        slong bits = FLINT_BITS - (slong) mod.norm;
        slong log_len = FLINT_BIT_COUNT(len2);

        if (2 * bits + log_len <= FLINT_BITS)
        {
            /* Set res[i] = poly1[i]*poly2[0] */
            mpn_mul_1(res, poly1, FLINT_MIN(len1, trunc), poly2[0]);

            if (len2 != 1)
            {
                /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
                if (trunc > len1)
                    mpn_mul_1(res + len1, poly2 + 1, trunc - len1,
                              poly1[len1 - 1]);

                /* out[i+j] += in1[i]*in2[j] */
                for (i = 0; i < FLINT_MIN(len1, trunc) - 1; i++)
                {
                    FLINT_ASSERT(FLINT_MIN(len2, trunc - i) > 1);
                    mpn_addmul_1(res + i + 1, poly2 + 1,
                                 FLINT_MIN(len2, trunc - i) - 1, poly1[i]);
                }
            }

            _nmod_vec_reduce(res, res, trunc, mod);
        }
        else
        {
            /* Set res[i] = poly1[i]*poly2[0] */
            _nmod_vec_scalar_mul_nmod(res, poly1, FLINT_MIN(len1, trunc),
                                 poly2[0], mod);

            if (len2 == 1)
                return;

            /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
            if (trunc > len1)
                _nmod_vec_scalar_mul_nmod(res + len1, poly2 + 1, trunc - len1,
                                     poly1[len1 - 1], mod);

            /* out[i+j] += in1[i]*in2[j] */
            for (i = 0; i < FLINT_MIN(len1, trunc) - 1; i++)
                _nmod_vec_scalar_addmul_nmod(res + i + 1, poly2 + 1,
                                        FLINT_MIN(len2, trunc - i) - 1, 
                                        poly1[i], mod);
        }
    }
}

void
nmod_poly_mullow_classical(nmod_poly_t res,
                           const nmod_poly_t poly1, const nmod_poly_t poly2,
                           slong trunc)
{
    slong len_out;

    if (poly1->length == 0 || poly2->length == 0 || trunc == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (trunc > len_out)
        trunc = len_out;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, trunc);
        if (poly1->length >= poly2->length)
            _nmod_poly_mullow_classical(temp->coeffs, poly1->coeffs,
                                        poly1->length, poly2->coeffs,
                                        poly2->length, trunc, poly1->mod);
        else
            _nmod_poly_mullow_classical(temp->coeffs, poly2->coeffs,
                                        poly2->length, poly1->coeffs,
                                        poly1->length, trunc, poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, trunc);
        if (poly1->length >= poly2->length)
            _nmod_poly_mullow_classical(res->coeffs, poly1->coeffs,
                                        poly1->length, poly2->coeffs,
                                        poly2->length, trunc, poly1->mod);
        else
            _nmod_poly_mullow_classical(res->coeffs, poly2->coeffs,
                                        poly2->length, poly1->coeffs,
                                        poly1->length, trunc, poly1->mod);
    }

    res->length = trunc;
    _nmod_poly_normalise(res);
}
