/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

/* Assumes poly1 and poly2 are not length 0. */
void
_nmod_poly_mul_classical(nn_ptr res, nn_srcptr poly1,
                         slong len1, nn_srcptr poly2, slong len2, nmod_t mod)
{
    slong i, j, n1, n2;
    int squaring;
    ulong c;

    if (len1 == 1)
    {
        res[0] = nmod_mul(poly1[0], poly2[0], mod);
        return;
    }

    if (len2 == 1)
    {
        _nmod_vec_scalar_mul_nmod(res, poly1, len1, poly2[0], mod);
        return;
    }

    squaring = (poly1 == poly2 && len1 == len2);

    const dot_params_t params = _nmod_vec_dot_params(FLINT_MIN(len1, len2), mod);

    if (params.method <= _DOT1)
    {
        flint_mpn_zero(res, len1 + len2 - 1);

        if (squaring)
        {
            for (i = 0; i < len1; i++)
            {
                c = poly1[i];

                res[2 * i] += c * c;
                c *= 2;

                for (j = i + 1; j < len1; j++)
                    res[i + j] += poly1[j] * c;
            }
        }
        else
        {
            for (i = 0; i < len1; i++)
            {
                ulong c = poly1[i];

                for (j = 0; j < len2; j++)
                    res[i + j] += c * poly2[j];
            }
        }

        _nmod_vec_reduce(res, res, len1 + len2 - 1, mod);
        return;
    }

    if (len2 == 2)
    {
        _nmod_vec_scalar_mul_nmod(res, poly1, len1, poly2[0], mod);
        _nmod_vec_scalar_addmul_nmod(res + 1, poly1, len1 - 1, poly2[1], mod);
        res[len1 + len2 - 2] = nmod_mul(poly1[len1 - 1], poly2[len2 - 1], mod);
        return;
    }

    if (squaring)
    {
        for (i = 0; i < 2 * len1 - 1; i++)
        {
            n1 = FLINT_MAX(0, i - len1 + 1);
            n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            c = _nmod_vec_dot_rev(poly1 + n1, poly1 + i - n2, n2 - n1 + 1, mod, params);
            c = nmod_add(c, c, mod);

            if (i % 2 == 0 && i / 2 < len1)
                NMOD_ADDMUL(c, poly1[i / 2], poly1[i / 2], mod);

            res[i] = c;
        }
    }
    else
    {
        for (i = 0; i < len1 + len2 - 1; i++)
        {
            n1 = FLINT_MIN(len1 - 1, i);
            n2 = FLINT_MIN(len2 - 1, i);

            res[i] = _nmod_vec_dot_rev(poly1 + i - n2,
                                       poly2 + i - n1,
                                       n1 + n2 - i + 1, mod, params);
        }
    }
}

void
nmod_poly_mul_classical(nmod_poly_t res,
                        const nmod_poly_t poly1, const nmod_poly_t poly2)
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
