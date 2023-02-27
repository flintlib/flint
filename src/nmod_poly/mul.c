/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

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

void _nmod_poly_mul(mp_ptr res, mp_srcptr poly1, slong len1, 
                             mp_srcptr poly2, slong len2, nmod_t mod)
{
    slong bits, cutoff_len;

    if (len2 <= 5)
    {
        _nmod_poly_mul_classical(res, poly1, len1, poly2, len2, mod);
        return;
    }

    bits = FLINT_BITS - (slong) mod.norm;
    cutoff_len = FLINT_MIN(len1, 2 * len2);

    if (3 * cutoff_len < 2 * FLINT_MAX(bits, 10))
        _nmod_poly_mul_classical(res, poly1, len1, poly2, len2, mod);
    else if (cutoff_len * bits < 800)
        _nmod_poly_mul_KS(res, poly1, len1, poly2, len2, 0, mod);
    else if (cutoff_len * (bits + 1) * (bits + 1) < 100000)
        _nmod_poly_mul_KS2(res, poly1, len1, poly2, len2, mod);
    else
        _nmod_poly_mul_KS4(res, poly1, len1, poly2, len2, mod);
}

void nmod_poly_mul(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    slong len1, len2, len_out;
    
    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_zero(res);

        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;

        nmod_poly_init2(temp, poly1->mod.n, len_out);

        if (len1 >= len2)
            _nmod_poly_mul(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, poly1->mod);
        else
            _nmod_poly_mul(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, poly1->mod);
        
        nmod_poly_swap(temp, res);
        nmod_poly_clear(temp);
    } else
    {
        nmod_poly_fit_length(res, len_out);
        
        if (len1 >= len2)
            _nmod_poly_mul(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, poly1->mod);
        else
            _nmod_poly_mul(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}
