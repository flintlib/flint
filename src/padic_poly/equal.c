/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "padic_poly.h"

int padic_poly_equal(const padic_poly_t poly1, const padic_poly_t poly2)
{
    if (poly1 == poly2)
    {
        return 1;
    }

    if (poly1->length != poly2->length || poly1->val != poly2->val)
    {
        return 0;
    }

    return _fmpz_vec_equal(poly1->coeffs, poly2->coeffs, poly1->length);
}
