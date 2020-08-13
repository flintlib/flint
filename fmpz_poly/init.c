/*
    Copyright (C) 2008, 2009 William Hart

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
fmpz_poly_init(fmpz_poly_t poly)
{
    poly->coeffs = NULL;

    poly->alloc = 0;
    poly->length = 0;
}

void
fmpz_poly_init2(fmpz_poly_t poly, slong alloc)
{
    if (alloc != 0)                  /* allocate space for alloc small coeffs */
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    else
        poly->coeffs = NULL;

    poly->alloc = alloc;
    poly->length = 0;
}
