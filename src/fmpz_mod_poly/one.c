/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_one(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, 1, ctx);
    fmpz_one(poly->coeffs + 0);
    _fmpz_mod_poly_set_length(poly, fmpz_is_one(fmpz_mod_ctx_modulus(ctx)) ? 0 : 1);
}
