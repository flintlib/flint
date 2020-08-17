/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

slong
fmpq_mat_rref_classical(fmpq_mat_t B, const fmpq_mat_t A)
{
    slong m, n, i, j, pivot_row, pivot_col, rank;

    m = A->r;
    n = A->c;

    if (m == 0 || n == 0)
        return 0;

    if (A != B)
        fmpq_mat_set(B, A);

    rank = 0;
    pivot_row = 0;
    pivot_col = 0;

    while (pivot_row < m && pivot_col < n)
    {
        if (!fmpq_mat_pivot(NULL, B, pivot_row, pivot_col))
        {
            pivot_col++;
            continue;
        }

        rank++;

        /* Scale row to have 1 as leading entry */
        for (j = pivot_col + 1; j < n; j++)
        {
            fmpq_div(fmpq_mat_entry(B, pivot_row, j),
                     fmpq_mat_entry(B, pivot_row, j),
                     fmpq_mat_entry(B, pivot_row, pivot_col));
        }

        /* Eliminate rows above and below */
        for (i = 0; i < m; i++)
        {
            if (i == pivot_row ||
                fmpq_is_zero(fmpq_mat_entry(B, i, pivot_col)))
                continue;

            for (j = pivot_col + 1; j < n; j++)
                fmpq_submul(fmpq_mat_entry(B, i, j),
                            fmpq_mat_entry(B, pivot_row, j),
                            fmpq_mat_entry(B, i, pivot_col));
        }

        /* Clear pivot column */
        for (i = 0; i < m; i++)
            fmpq_set_si(fmpq_mat_entry(B, i, pivot_col), i == pivot_row, 1);

        pivot_row++;
        pivot_col++;
    }

    return rank;
}
