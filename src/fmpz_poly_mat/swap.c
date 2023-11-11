/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_swap(fmpz_poly_mat_t A, fmpz_poly_mat_t B)
{
    FLINT_SWAP(fmpz_poly_mat_struct, *A, *B);
}

void
fmpz_poly_mat_swap_entrywise(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < fmpz_poly_mat_nrows(mat1); i++)
        for (j = 0; j < fmpz_poly_mat_ncols(mat1); j++)
            fmpz_poly_swap(fmpz_poly_mat_entry(mat2, i, j), fmpz_poly_mat_entry(mat1, i, j));
}
