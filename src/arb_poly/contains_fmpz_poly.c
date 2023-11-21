/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

int
arb_poly_contains_fmpz_poly(const arb_poly_t poly1, const fmpz_poly_t poly2)
{
    slong i;

    if (poly2->length > poly1->length)
        return 0;

    for (i = 0; i < poly2->length; i++)
        if (!arb_contains_fmpz(poly1->coeffs + i, poly2->coeffs + i))
            return 0;

    for (i = poly2->length; i < poly1->length; i++)
        if (!arb_contains_zero(poly1->coeffs + i))
            return 0;

    return 1;
}

