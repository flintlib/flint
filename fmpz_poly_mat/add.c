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
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_add(fmpz_poly_mat_t C,
                        const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_poly_add(fmpz_poly_mat_entry(C, i, j),
                          fmpz_poly_mat_entry(A, i, j),
                          fmpz_poly_mat_entry(B, i, j));
}
