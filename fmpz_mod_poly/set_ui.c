/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_set_ui(fmpz_mod_poly_t f, ulong x, const fmpz_mod_ctx_t ctx)
{
    if (x == 0)
    {
        fmpz_mod_poly_zero(f, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(f, 1, ctx);
        _fmpz_mod_poly_set_length(f, 1);
        fmpz_set_ui(f->coeffs, x);
        fmpz_mod(f->coeffs, f->coeffs, fmpz_mod_ctx_modulus(ctx));
        _fmpz_mod_poly_normalise(f);
    }
}
