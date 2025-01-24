/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

int
fmpq_mat_pivot(slong * perm, fmpq_mat_t mat, slong r, slong c)
{
    slong j;

    if (!fmpq_is_zero(fmpq_mat_entry(mat, r, c)))
        return 1;

    for (j = r + 1; j < mat->r; j++)
    {
        if (!fmpq_is_zero(fmpq_mat_entry(mat, j, c)))
        {
            fmpq_mat_swap_rows(mat, perm, j, r);
            return -1;
        }
    }
    return 0;
}
