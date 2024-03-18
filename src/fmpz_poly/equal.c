/*
    Copyright (C) 2008, 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

int
fmpz_poly_equal(const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    slong i;

    if (poly1 == poly2)
        return 1;               /* same polynomial */

    if (poly1->length != poly2->length)
        return 0;               /* check if lengths the same */

    for (i = 0; i < poly1->length; i++) /* check if coefficients the same */
        if (!fmpz_equal(poly1->coeffs + i, poly2->coeffs + i))
            return 0;

    return 1;
}

int
fmpz_poly_equal_fmpz(const fmpz_poly_t poly, const fmpz_t c)
{
	return ((poly->length == 0) && fmpz_is_zero(c)) ||
        ((poly->length == 1) && fmpz_equal(poly->coeffs, c));
}
