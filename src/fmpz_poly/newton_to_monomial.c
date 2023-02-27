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
_fmpz_poly_newton_to_monomial(fmpz * poly, const fmpz * roots, slong n)
{
    slong i, j;

    for (i = n - 2; i >= 0; i--)
        for (j = i; j < n - 1; j++)
            fmpz_submul(poly + j, poly + j + 1, roots + i);
}
