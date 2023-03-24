/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_mat.h"

void
fmpz_mod_mat_swap(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2)
{
    if (mat1 != mat2)
    {
        fmpz_mod_mat_struct tmp;

        tmp = *mat1;
        *mat1 = *mat2;
        *mat2 = tmp;
    }
}

void
fmpz_mod_mat_swap_entrywise(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < fmpz_mod_mat_nrows(mat1); i++)
        for (j = 0; j < fmpz_mod_mat_ncols(mat1); j++)
            fmpz_swap(fmpz_mod_mat_entry(mat2, i, j), fmpz_mod_mat_entry(mat1, i, j));
}
