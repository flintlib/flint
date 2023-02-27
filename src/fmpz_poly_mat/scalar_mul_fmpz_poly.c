/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_scalar_mul_fmpz_poly(fmpz_poly_mat_t B, const fmpz_poly_mat_t A,
                                                        const fmpz_poly_t c)
{
    slong i, j;

    for (i = 0; i < fmpz_poly_mat_nrows(B); i++)
        for (j = 0; j < fmpz_poly_mat_ncols(B); j++)
            fmpz_poly_mul(fmpz_poly_mat_entry(B, i, j),
                          fmpz_poly_mat_entry(A, i, j), c);
}
