/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mini.h"
#include "fmpz_mat.h"

void
fmpz_mat_swap_cols(fmpz_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s && !fmpz_mat_is_empty(mat))
    {
        slong t;

        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

       for (t = 0; t < mat->r; t++)
       {
           fmpz_swap(fmpz_mat_entry(mat, t, r), fmpz_mat_entry(mat, t, s));
       }
    }
}

