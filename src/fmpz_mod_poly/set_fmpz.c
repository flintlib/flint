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
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_set_fmpz(fmpz_mod_poly_t poly, const fmpz_t c,
                                                      const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(poly, 1, ctx);
    fmpz_mod(poly->coeffs, c, fmpz_mod_ctx_modulus(ctx));
    _fmpz_mod_poly_set_length(poly, 1);
    _fmpz_mod_poly_normalise(poly);
}

