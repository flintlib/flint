/*
    Copyright (C) 2011 Sebastian Pancratz

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
#include "fmpz_mod_poly.h"


void fmpz_mod_poly_init2(fmpz_mod_poly_t poly, slong alloc,
                                                      const fmpz_mod_ctx_t ctx)
{
    if (alloc)                  /* allocate space for alloc small coeffs */
        poly->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
    else
        poly->coeffs = NULL;

    poly->alloc = alloc;
    poly->length = 0;
}

