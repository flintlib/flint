/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "padic_mat.h"

void padic_mat_swap(padic_mat_t A, padic_mat_t B)
{
    FLINT_SWAP(padic_mat_struct, *A, *B);
}

void padic_mat_swap_entrywise(padic_mat_t mat1, padic_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < padic_mat_nrows(mat1); i++)
        for (j = 0; j < padic_mat_ncols(mat1); j++)
            fmpz_swap(padic_mat_entry(mat2, i, j), padic_mat_entry(mat1, i, j));
}
