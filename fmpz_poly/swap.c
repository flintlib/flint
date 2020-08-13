/*
    Copyright (C) 2008, 2009 William Hart

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
fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2)
{
    if (poly1 != poly2)
    {
        slong temp;
        fmpz *temp_c;

        temp = poly1->length;
        poly1->length = poly2->length;
        poly2->length = temp;

        temp = poly1->alloc;
        poly1->alloc = poly2->alloc;
        poly2->alloc = temp;

        temp_c = poly1->coeffs;
        poly1->coeffs = poly2->coeffs;
        poly2->coeffs = temp_c;
    }
}
