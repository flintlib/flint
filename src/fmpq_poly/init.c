/*
    Copyright (C) 2010 William Hart
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
#include "fmpq_poly.h"

void fmpq_poly_init(fmpq_poly_t poly)
{
    poly->coeffs = NULL;
    fmpz_init(poly->den);
    fmpz_one(poly->den);
    poly->alloc  = 0;
    poly->length = 0;
}

void fmpq_poly_init2(fmpq_poly_t poly, slong alloc)
{
    /* Allocate space for alloc small coeffs */
    poly->coeffs = (alloc ? (fmpz *) flint_calloc(alloc, sizeof(fmpz)) : NULL);
    
    fmpz_init(poly->den);
    fmpz_one(poly->den);
    poly->alloc  = alloc;
    poly->length = 0;
}

