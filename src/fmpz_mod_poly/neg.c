/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_neg(fmpz *res, const fmpz *poly, slong len, const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_vec_neg(res, poly, len, ctx);
}

void fmpz_mod_poly_neg(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong len = poly->length;

    fmpz_mod_poly_fit_length(res, len, ctx);
    _fmpz_mod_poly_set_length(res, len);

    _fmpz_mod_poly_neg(res->coeffs, poly->coeffs, poly->length, ctx);
}
