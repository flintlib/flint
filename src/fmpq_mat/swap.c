/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

void
fmpq_mat_swap(fmpq_mat_t mat1, fmpq_mat_t mat2)
{
    FLINT_SWAP(fmpq_mat_struct, *mat1, *mat2);
}

void
fmpq_mat_swap_entrywise(fmpq_mat_t mat1, fmpq_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < fmpq_mat_nrows(mat1); i++)
        for (j = 0; j < fmpq_mat_ncols(mat1); j++)
            fmpq_swap(fmpq_mat_entry(mat2, i, j), fmpq_mat_entry(mat1, i, j));
}

void
fmpq_mat_swap_cols(fmpq_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s && !fmpq_mat_is_empty(mat))
    {
        slong i;

        if (perm)
            FLINT_SWAP(slong, perm[r], perm[s]);

        for (i = 0; i < mat->r; i++)
            fmpq_swap(fmpq_mat_entry(mat, i, r), fmpq_mat_entry(mat, i, s));
    }
}

void
fmpq_mat_swap_rows(fmpq_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s && !fmpq_mat_is_empty(mat))
    {
        if (perm)
            FLINT_SWAP(slong, perm[r], perm[s]);

        FLINT_SWAP(fmpq *, mat->rows[r], mat->rows[s]);
    }
}
