/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_swap(fmpz_mat_t mat1, fmpz_mat_t mat2)
{
    FLINT_SWAP(fmpz_mat_struct, *mat1, *mat2);
}

void
fmpz_mat_swap_entrywise(fmpz_mat_t mat1, fmpz_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < fmpz_mat_nrows(mat1); i++)
        for (j = 0; j < fmpz_mat_ncols(mat1); j++)
            fmpz_swap(fmpz_mat_entry(mat2, i, j), fmpz_mat_entry(mat1, i, j));
}

void
fmpz_mat_swap_cols(fmpz_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s && !fmpz_mat_is_empty(mat))
    {
        slong i;

        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

       for (i = 0; i < mat->r; i++)
           fmpz_swap(fmpz_mat_entry(mat, i, r), fmpz_mat_entry(mat, i, s));
    }
}
