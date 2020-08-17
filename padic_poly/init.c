/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011, 2012 Sebastian Pancratz

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
#include "padic_poly.h"

void padic_poly_init(padic_poly_t poly)
{
    poly->coeffs = NULL;
    poly->alloc  = 0;
    poly->length = 0;
    poly->val    = 0;
    poly->N      = PADIC_DEFAULT_PREC;
}

void padic_poly_init2(padic_poly_t poly, slong alloc, slong prec)
{
    poly->coeffs = alloc ? (fmpz *) flint_calloc(alloc, sizeof(fmpz)) : NULL;
    poly->alloc  = alloc;
    poly->length = 0;
    poly->val    = 0;
    poly->N      = prec;
}

