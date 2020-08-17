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
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"

void fmpq_poly_swap(fmpq_poly_t poly1, fmpq_poly_t poly2)
{
    slong t;
    fmpz * tptr;
    
    t             = poly1->length;
    poly1->length = poly2->length;
    poly2->length = t;
    
    t             = poly1->alloc;
    poly1->alloc  = poly2->alloc;
    poly2->alloc  = t;
    
    tptr          = poly1->coeffs;
    poly1->coeffs = poly2->coeffs;
    poly2->coeffs = tptr;
    
    fmpz_swap(poly1->den, poly2->den);
}

