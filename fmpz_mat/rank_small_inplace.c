/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "perm.h"
#include "longlong.h"

#define E(j,k) fmpz_mat_entry(B,j,k)

slong _fmpz_mat_rank_overflow(fmpz_mat_t B, slong pivot_row, slong pivot_col)
{
    fmpz_t den;
    fmpz_mat_t window;
    slong rank, m = B->r, n = B->c;

    fmpz_init(den);
    fmpz_mat_window_init(window, B, pivot_row, pivot_col, m, n);

    if (FLINT_MIN(m - pivot_row, n - pivot_col) < 25)
        rank = fmpz_mat_fflu(window, den, NULL, window, 0);
    else
        rank = fmpz_mat_rref(window, den, window);

    fmpz_mat_window_clear(window);
    fmpz_clear(den);

    return rank;
}

slong
fmpz_mat_rank_small_inplace(fmpz_mat_t B)
{
    slong m, n, j, k, rank, q, r = -1, pivot_row, pivot_col;
    /* mask to check that a*b + c will be small */
    const ulong mask = ~((UWORD(1)<<((FLINT_BITS - 4)/2)) - 1);
    ulong largest;

    if (fmpz_mat_is_empty(B))
    {
        return 0;
    }

    m = B->r;
    n = B->c;
    rank = pivot_row = pivot_col = 0;

    if (pivot_row < m && pivot_col < n)
       r = fmpz_mat_find_pivot_smallest(B, pivot_row, m, pivot_col);

    while (pivot_row < m && pivot_col < n)
    {
        if (r == -1)
        {
            pivot_col++;
            if (pivot_col == n)
                break;

            r = fmpz_mat_find_pivot_smallest(B, pivot_row, m, pivot_col);

            continue;
        }
        else if (r != pivot_row)
        {
            fmpz_mat_swap_rows(B, NULL, pivot_row, r);
        }

        largest = 0;

        for (j = pivot_row + 1; j < m; j++)
        {
            if (*E(j, pivot_col) != 0)
            {
                q = (*E(j, pivot_col))/(*E(pivot_row, pivot_col));

                for (k = pivot_col; k < n; k++)
                {
                    (*E(j, k)) -= q*(*E(pivot_row, k));
                    largest |= FLINT_ABS(*E(j, k));
                }
            }
        }

        if ((largest & mask) != 0) /* rest may overflow, call safe code */
            return rank + _fmpz_mat_rank_overflow(B, pivot_row, pivot_col);

        r = fmpz_mat_find_pivot_smallest(B, pivot_row + 1, m, pivot_col);
        
        if (r == -1)
        {
            pivot_row++;
            pivot_col++;

            rank++;

            if (pivot_row < m && pivot_col < n)
               r = fmpz_mat_find_pivot_smallest(B, pivot_row, m, pivot_col);
        }
    }

    return rank;
}
