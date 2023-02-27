/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
fmpz_poly_clear(fmpz_poly_t poly)
{
    if (poly->coeffs)
    {
        slong i;
        for (i = 0; i < poly->alloc; i++)
            _fmpz_demote(poly->coeffs + i);
        flint_free(poly->coeffs);
    }
}
