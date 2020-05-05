/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

slong
fmpz_poly_remove(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    fmpz_poly_t q, p;
    slong i = 0;

    if (poly2->length == 0)
    {
        flint_printf("Exception (fmpz_poly_remove). Division by zero.\n");
	flint_abort();
    }

    if (poly2->length == 1 && fmpz_is_pm1(poly2->coeffs + 0))
    {
        flint_printf("Exception (fmpz_poly_remove). Divisor must not be a unit.\n");
	flint_abort();
    }

    if (poly2->length > poly1->length)
    {
        fmpz_poly_set(res, poly1);
	return 0;
    }

    fmpz_poly_init(q);
    fmpz_poly_init(p);

    fmpz_poly_set(p, poly1);
    
    while (fmpz_poly_divides(q, p, poly2))
    {
        fmpz_poly_swap(p, q);
	i++;
    }

    fmpz_poly_set(res, p);

    fmpz_poly_clear(p);
    fmpz_poly_clear(q);

    return i;
}
