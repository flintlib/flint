/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_invert_cols(fmpz_mat_t mat, slong * perm)
{
    if (!fmpz_mat_is_empty(mat))
    {
        slong t;
        slong i;
        slong c = mat->c;
        slong k = mat->c/2;

        if (perm)
        {
            for (i =0; i < k; i++)
            {
                t = perm[i];
                perm[i] = perm[c - i];
                perm[c - i] = t;
            }
        }

        for (t = 0; t < mat->r; t++)
        {
            for (i = 0; i < k; i++)
            {
                fmpz_swap(fmpz_mat_entry(mat, t, i), fmpz_mat_entry(mat, t, c - i - 1));
            }
        }
    }
}

