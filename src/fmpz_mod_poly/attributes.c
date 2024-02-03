/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

/* bogus for non-prime modulus */
int fmpz_mod_poly_is_unit(const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx)
{
    return (op->length == 1) && fmpz_mod_is_invertible(op->coeffs + 0, ctx);
}
