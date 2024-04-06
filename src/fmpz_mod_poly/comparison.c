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

int fmpz_mod_poly_equal(const fmpz_mod_poly_t poly1,
                        const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
{
    return fmpz_poly_equal((fmpz_poly_struct *) poly1,
                           (fmpz_poly_struct *) poly2);
}

int fmpz_mod_poly_equal_trunc(const fmpz_mod_poly_t poly1,
                const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
{
    return fmpz_poly_equal_trunc((fmpz_poly_struct *) poly1,
                           (fmpz_poly_struct *) poly2, n);
}

int fmpz_mod_poly_is_zero(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t FLINT_UNUSED(ctx))
{
    return poly->length == 0;
}
