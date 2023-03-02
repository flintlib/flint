/*
    Copyright (C) 2011 Fredrik Johansson

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
fmpz_poly_evaluate_fmpz_vec(fmpz * res, const fmpz_poly_t f,
                                const fmpz * a, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        fmpz_poly_evaluate_fmpz(res + i, f, a + i);
}
