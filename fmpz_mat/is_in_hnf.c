/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

int fmpz_mat_is_in_hnf(const fmpz_mat_t A)
{
    slong i, last_i, j, prev_j;

    /* find last non-zero row */
    for (last_i = A->r - 1; last_i >= 0; last_i--)
    {
        for (j = 0; j < A->c; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, last_i, j)))
                break;
        }
        if (j < A->c)
            break;
    }

    /* hermite form structure */
    prev_j = -1;
    for (i = 0; i <= last_i; i++)
    {
        slong i2;

        for (j = 0; j < A->c; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
            {
                if (fmpz_sgn(fmpz_mat_entry(A, i, j)) < 0)
                    return 0;
                break;
            }
        }
        if (j == A->c || j <= prev_j)
            return 0;
        prev_j = j;
        for (i2 = 0; i2 < i; i2++)
        {
            if (fmpz_cmp(fmpz_mat_entry(A, i2, j), fmpz_mat_entry(A, i, j)) >= 0)
                return 0;
            if (fmpz_sgn(fmpz_mat_entry(A, i2, j)) < 0)
                return 0;
        }
    }

    return 1;
}
