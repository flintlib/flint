/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

char *
fmpz_mod_poly_get_str(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
{
    return _fmpz_poly_get_str(poly->coeffs, poly->length);
}

char *
fmpz_mod_poly_get_str_pretty(const fmpz_mod_poly_t poly, const char * x,
                                                      const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
{
    return _fmpz_poly_get_str_pretty(poly->coeffs, poly->length, x);
}
