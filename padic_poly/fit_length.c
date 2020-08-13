/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>

#include "fmpz.h"
#include "padic_poly.h"

void padic_poly_fit_length(padic_poly_t poly, slong len)
{
    if (len > poly->alloc)
    {
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;

        if (poly->alloc)           /* Realloc */
        {
            poly->coeffs = (fmpz *) flint_realloc(poly->coeffs, len * sizeof(fmpz));
            mpn_zero((mp_ptr) (poly->coeffs + poly->alloc), len - poly->alloc);
        }
        else                       /* Nothing allocated already so do it now */
        {
            poly->coeffs = (fmpz *) flint_calloc(len, sizeof(fmpz));
        }

        poly->alloc = len;
    }
}

