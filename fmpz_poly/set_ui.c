/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
fmpz_poly_set_ui(fmpz_poly_t poly, ulong c)
{
    if (c == UWORD(0))
        fmpz_poly_zero(poly);
    else
    {
        fmpz_poly_fit_length(poly, 1);
        fmpz_set_ui(poly->coeffs, c);
        _fmpz_poly_set_length(poly, 1);
    }
}
