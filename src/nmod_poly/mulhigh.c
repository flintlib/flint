/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_mulhigh(mp_ptr res, mp_srcptr poly1, slong len1,
                             mp_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    slong bits, bits2;

    if (len1 + len2 <= 6)
    {
        _nmod_poly_mulhigh_classical(res, poly1, len1, poly2, len2, n, mod);
        return;
    }

    bits = FLINT_BITS - (slong) mod.norm;
    bits2 = FLINT_BIT_COUNT(len1);

    if (2 * bits + bits2 <= FLINT_BITS && len1 + len2 < 16)
        _nmod_poly_mulhigh_classical(res, poly1, len1, poly2, len2, n, mod);
    else
        _nmod_poly_mul_KS(res, poly1, len1, poly2, len2, 0, mod);
}

void nmod_poly_mulhigh(nmod_poly_t res,
                   const nmod_poly_t poly1, const nmod_poly_t poly2, slong start)
{
    slong len1, len2, len_out;

    len1 = poly1->length;
    len2 = poly2->length;
    len_out = len1 + len2 - 1;

    if (len1 == 0 || len2 == 0 || start >= len_out)
    {
        nmod_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;

        nmod_poly_init2(temp, poly1->mod.n, len_out);

        if (len1 >= len2)
            _nmod_poly_mulhigh(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, start, poly1->mod);
        else
            _nmod_poly_mulhigh(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, start, poly1->mod);

        nmod_poly_swap(temp, res);
        nmod_poly_clear(temp);
    } else
    {
        nmod_poly_fit_length(res, len_out);

        if (len1 >= len2)
            _nmod_poly_mulhigh(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, start, poly1->mod);
        else
            _nmod_poly_mulhigh(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, start, poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}

/* Assumes poly1 and poly2 are not length 0. */
void
_nmod_poly_mulhigh_classical(mp_ptr res, mp_srcptr poly1,
                             slong len1, mp_srcptr poly2, slong len2, slong start,
                             nmod_t mod)
{
    slong m, n;

    _nmod_vec_zero(res, start);

    if (len1 == 1)              /* Special case if the length of both inputs is 1 */
    {
        if (start == 0)
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
            if (start < len1)
                mpn_mul_1(res + start, poly1 + start, len1 - start, poly2[0]);

            if (len2 != 1)
            {
                /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
                m = FLINT_MAX(len1 - 1, start);

                FLINT_ASSERT(len2 - 1 + len1 - m > 0);
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
                            slong start)
{
    slong len_out = poly1->length + poly2->length - 1;

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
