/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

void fmpz_mod_poly_mulhigh(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
            const fmpz_mod_poly_t poly2, slong start, const fmpz_mod_ctx_t ctx)
{
    slong len1, len2, len_out;
    
    len1 = poly1->length;
    len2 = poly2->length;
    len_out = len1 + len2 - 1;

    if (start == 0)
    {
        fmpz_mod_poly_mul(res, poly1, poly2, ctx);
        return;
    }

    if (len1 < 1 || len2 < 1 || start >= len_out)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_mod_poly_t temp;
        fmpz_mod_poly_init2(temp, len_out, ctx);

        if (len1 >= len2)
            _fmpz_poly_mulhigh(temp->coeffs, poly1->coeffs, len1,
                                         poly2->coeffs, len2, start);
        else
            _fmpz_poly_mulhigh(temp->coeffs, poly2->coeffs, len2,
                                         poly1->coeffs, len1, start);
        
        fmpz_mod_poly_swap(temp, res, ctx);
        fmpz_mod_poly_clear(temp, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, len_out, ctx);
        
        if (len1 >= len2)
            _fmpz_poly_mulhigh(res->coeffs, poly1->coeffs, len1,
                                         poly2->coeffs, len2, start);
        else
            _fmpz_poly_mulhigh(res->coeffs, poly2->coeffs, len2,
                                         poly1->coeffs, len1, start);
    }

    _fmpz_vec_scalar_mod_fmpz(res->coeffs, res->coeffs, len_out,
                                                    fmpz_mod_ctx_modulus(ctx));
    _fmpz_mod_poly_set_length(res, len_out);
    _fmpz_mod_poly_normalise(res);
}
