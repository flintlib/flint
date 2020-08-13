/*
    Copyright (C) 2010-2012 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

int fmpz_mat_is_in_rref_with_rank(const fmpz_mat_t A, const fmpz_t den, slong rank)
{
    slong i, j, k, prev_pivot;

    /* bottom should be zero */
    for (i = rank; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
                return 0;

    prev_pivot = -1;

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
            {
                /* pivot should have a higher column index than previous */
                if (j <= prev_pivot)
                    return 0;

                /* column should be 0 ... 0 1 0 ... 0 */
                for (k = 0; k < rank; k++)
                {
                    if (i == k && !fmpz_equal(fmpz_mat_entry(A, k, j), den))
                        return 0;
                    if (i != k && !fmpz_is_zero(fmpz_mat_entry(A, k, j)))
                        return 0;
                }

                prev_pivot = j;
                break;
            }
        }
    }

    return 1;
}
