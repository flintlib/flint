/*
    Copyright (C) 2011, 2010 Sebastian Pancratz

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

void fmpz_mod_poly_zero_coeffs(fmpz_mod_poly_t poly, slong i, slong j,
                                                      const fmpz_mod_ctx_t ctx)
{
    if (i < 0)
        i = 0;
    if (j > poly->length)
        j = poly->length;

    _fmpz_vec_zero(poly->coeffs + i, j - i);

    if (j == poly->length)
    {
        _fmpz_mod_poly_set_length(poly, i);
        _fmpz_mod_poly_normalise(poly);
    }
}

