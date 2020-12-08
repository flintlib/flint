/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_gen(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, 2, ctx);
    fmpz_zero(poly->coeffs + 0);
    fmpz_one(poly->coeffs + 1);
    _fmpz_mod_poly_set_length(poly, 2*!fmpz_is_one(fmpz_mod_ctx_modulus(ctx)));
}

