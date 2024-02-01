/*
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
fmpz_poly_fit_length(fmpz_poly_t poly, slong len)
{
    if (len > poly->alloc)
    {
        /* At least double number of allocated coeffs */
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;
        fmpz_poly_realloc(poly, len);
    }
}
