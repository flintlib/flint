/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void _fmpz_poly_lcm(fmpz * res, const fmpz * poly1, slong len1, 
                                              const fmpz * poly2, slong len2)
{
    fmpz *W;
    slong lenW = len2;

    W = _fmpz_vec_init(len2);

    _fmpz_poly_mul(res, poly1, len1, poly2, len2);

    _fmpz_poly_gcd(W, poly1, len1, poly2, len2);

    FMPZ_VEC_NORM(W, lenW);

    if (lenW == 1)
    {
        if (fmpz_sgn(res + (len1 + len2 - 1 - 1)) < 0)
            fmpz_neg(W + 0, W + 0);
        _fmpz_vec_scalar_divexact_fmpz(res, res, len1 + len2 - 1, W + 0);
    }
    else
    {
        fmpz *V;
        slong lenV = len1 + len2 - lenW;

        V = _fmpz_vec_init(lenV);
        _fmpz_poly_div(V, res, len1 + len2 - 1, W, lenW, 0);
        if (fmpz_sgn(V + (lenV - 1)) > 0)
            _fmpz_vec_set(res, V, lenV);
        else
            _fmpz_vec_neg(res, V, lenV);
        _fmpz_vec_zero(res + lenV, len1 + len2 - 1 - lenV);
        _fmpz_vec_clear(V, lenV);
    }

    _fmpz_vec_clear(W, len2);
}

void fmpz_poly_lcm(fmpz_poly_t res, const fmpz_poly_t poly1, 
                                                    const fmpz_poly_t poly2)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;

        fmpz_poly_init(t);
        fmpz_poly_lcm(t, poly1, poly2);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    fmpz_poly_fit_length(res, len1 + len2 - 1);
    _fmpz_poly_set_length(res, len1 + len2 - 1);
    if (len1 >= len2)
    {
        _fmpz_poly_lcm(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2);
    }
    else
    {
        _fmpz_poly_lcm(res->coeffs, poly2->coeffs, len2, poly1->coeffs, len1);
    }
    _fmpz_poly_normalise(res);
    
}

