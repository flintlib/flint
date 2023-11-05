/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

void fmpq_mat_invert_rows(fmpq_mat_t mat, slong * perm)
{
    slong i;

    for (i = 0; i < mat->r/2; i++)
        fmpq_mat_swap_rows(mat, perm, i, mat->r - i - 1);
}

void fmpq_mat_invert_cols(fmpq_mat_t mat, slong * perm)
{
    if (!fmpq_mat_is_empty(mat))
    {
        slong t, i;
        slong c = mat->c;
        slong k = mat->c/2;

        if (perm != NULL)
            for (i = 0; i < k; i++)
                FLINT_SWAP(slong, perm[i], perm[c - i - 1]);

        for (t = 0; t < mat->r; t++)
            for (i = 0; i < k; i++)
                fmpq_swap(fmpq_mat_entry(mat, t, i), fmpq_mat_entry(mat, t, c - i - 1));
    }
}
