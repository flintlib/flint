/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

#define E(j,k) fmpz_mat_entry(A,j,k)

len_t
fmpz_mat_rref_mod(len_t *perm, fmpz_mat_t A, const fmpz_t p)
{
    fmpz_t t, inv;
    len_t m, n, j, k, rank, r, pivot_row, pivot_col;

    if (fmpz_mat_is_empty(A))
    {
        return 0;
    }

    m = A->r;
    n = A->c;
    rank = pivot_row = pivot_col = 0;

    fmpz_init(t);
    fmpz_init(inv);

    while (pivot_row < m && pivot_col < n)
    {
        r = fmpz_mat_find_pivot_any(A, pivot_row, m, pivot_col);

        if (r == -1)
        {
            pivot_col++;
            continue;
        }
        else if (r != pivot_row)
        {
            fmpz_mat_swap_rows(A, perm, pivot_row, r);
        }
        rank++;

        fmpz_invmod(inv, E(pivot_row, pivot_col), p);

        /* pivot row */
        for (k = pivot_col + 1; k < n; k++)
        {
            fmpz_mul(E(pivot_row, k), E(pivot_row, k), inv);
            fmpz_mod(E(pivot_row, k), E(pivot_row, k), p);
        }
        fmpz_one(E(pivot_row, pivot_col));

        /* other rows */
        for (j = 0; j < m; j++)
        {
            if (j == pivot_row)
                continue;

            for (k = pivot_col + 1; k < n; k++)
            {
                fmpz_mul(t, E(j, pivot_col), E(pivot_row, k));
                fmpz_sub(E(j, k), E(j, k), t);
                fmpz_mod(E(j, k), E(j, k), p);
            }

            fmpz_zero(E(j, pivot_col));
        }

        pivot_row++;
        pivot_col++;
    }

    fmpz_clear(inv);
    fmpz_clear(t);

    return rank;
}
