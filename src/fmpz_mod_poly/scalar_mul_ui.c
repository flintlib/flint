/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_scalar_mul_ui(fmpz *res, const fmpz *poly, slong len,
                                    ulong x, const fmpz_mod_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mod_set_ui(t, x, ctx);
    _fmpz_mod_vec_scalar_mul_fmpz_mod(res, poly, len, t, ctx);
    fmpz_clear(t);
}

void fmpz_mod_poly_scalar_mul_ui(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, ulong x, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(res, poly->length, ctx);
    _fmpz_mod_poly_scalar_mul_ui(res->coeffs, poly->coeffs, poly->length, x, ctx);
    _fmpz_mod_poly_set_length(res, poly->length);
    _fmpz_mod_poly_normalise(res);
}

