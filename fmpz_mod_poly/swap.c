/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_swap(fmpz_mod_poly_t poly1, fmpz_mod_poly_t poly2)
{
    if (poly1 != poly2)
    {
        slong t;
        fmpz *c;

        t             = poly1->length;
        poly1->length = poly2->length;
        poly2->length = t;

        t            = poly1->alloc;
        poly1->alloc = poly2->alloc;
        poly2->alloc = t;

        c             = poly1->coeffs;
        poly1->coeffs = poly2->coeffs;
        poly2->coeffs = c;
    }
}

