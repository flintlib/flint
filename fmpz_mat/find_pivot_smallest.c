/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

slong
fmpz_mat_find_pivot_smallest(const fmpz_mat_t mat,
                                    slong start_row, slong end_row, slong c)
{
    slong r, smallest = -1;
    fmpz * currptr, * smallptr = NULL;

    for (r = start_row; r < end_row; r++)
    {
        currptr = fmpz_mat_entry(mat, r, c);
        if (!fmpz_is_zero(currptr))
        {
            if (0 > smallest || fmpz_cmpabs(currptr, smallptr) < 0)
            {
                smallest = r;
                smallptr = currptr;
            }
        }
    }

    return smallest;
}
