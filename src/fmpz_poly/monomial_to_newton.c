/*
    Copyright (C) 2012 Fredrik Johansson

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
_fmpz_poly_monomial_to_newton(fmpz * poly, const fmpz * roots, slong n)
{
    slong i, j;

    for (i = 0; i < n - 1; i++)
        for (j = n - 2; j >= i; j--)
            fmpz_addmul(poly + j, poly + j + 1, roots + i);
}
