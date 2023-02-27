/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

void padic_poly_swap(padic_poly_t poly1, padic_poly_t poly2)
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

        t          = poly1->val;
        poly1->val = poly2->val;
        poly2->val = t;

        t        = poly1->N;
        poly1->N = poly2->N;
        poly2->N = t;

        c             = poly1->coeffs;
        poly1->coeffs = poly2->coeffs;
        poly2->coeffs = c;
    }
}

