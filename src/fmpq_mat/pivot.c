/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

int
fmpq_mat_pivot(slong * perm, fmpq_mat_t mat, slong r, slong c)
{
    slong t, j;
    fmpq * u;

    if (!fmpq_is_zero(fmpq_mat_entry(mat, r, c)))
        return 1;

    for (j = r + 1; j < mat->r; j++)
    {
        if (!fmpq_is_zero(fmpq_mat_entry(mat, j, c)))
        {
            if (perm)
            {
                t = perm[j];
                perm[j] = perm[r];
                perm[r] = t;
            }

            u = mat->rows[j];
            mat->rows[j] = mat->rows[r];
            mat->rows[r] = u; 
            return -1;
        }
    }
    return 0;
}
