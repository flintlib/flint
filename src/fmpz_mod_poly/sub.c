/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_sub(fmpz *res, const fmpz *poly1, slong len1,
                                   const fmpz *poly2, slong len2, const fmpz_mod_ctx_t ctx)
{
    slong len3 = FLINT_MIN(len1, len2);

    _fmpz_mod_vec_sub(res, poly1, poly2, len3, ctx);

    if (len1 > len3)
        _fmpz_vec_set(res + len3, poly1 + len3, len1 - len3);
    if (len2 > len3)
        _fmpz_mod_vec_neg(res + len3, poly2 + len3, len2 - len3, ctx);
}

void fmpz_mod_poly_sub(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                         const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    fmpz_mod_poly_fit_length(res, max, ctx);
    _fmpz_mod_poly_sub(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, ctx);
    _fmpz_mod_poly_set_length(res, max);
    _fmpz_mod_poly_normalise(res);
}
