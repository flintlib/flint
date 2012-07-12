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

******************************************************************************/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

#define E(j,k) fmpz_mat_entry(B,j,k)

long
fmpz_mat_rref(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
{
    fmpz_t t;
    long m, n, j, k, rank, r, pivot_row, pivot_col;
    int sign;

    if (fmpz_mat_is_empty(A))
    {
        fmpz_set_ui(den, 1);
        return 0;
    }

    fmpz_mat_set(B, A);
    m = B->r;
    n = B->c;
    rank = pivot_row = pivot_col = 0;
    sign = 1;

    fmpz_init(t);

    while (pivot_row < m && pivot_col < n)
    {
        r = fmpz_mat_find_pivot_any(B, pivot_row, m, pivot_col);

        if (r == -1)
        {
            pivot_col++;
            continue;
        }
        else if (r != pivot_row)
        {
            fmpz_mat_swap_rows(B, NULL, pivot_row, r);
            sign = -sign;
        }
        rank++;

        for (j = 0; j < m; j++)
        {
            if (j == pivot_row)
                continue;

            for (k = j > pivot_row ? pivot_col + 1 : j; k < n; k++)
            {
                if (k == pivot_col)
                    continue;

                fmpz_mul(E(j, k), E(j, k), E(pivot_row, pivot_col));
                fmpz_mul(t, E(j, pivot_col), E(pivot_row, k));
                fmpz_sub(E(j, k), E(j, k), t);

                if (pivot_row > 0)
                    fmpz_divexact(E(j, k), E(j, k), den);
            }

            fmpz_zero(E(j, pivot_col));
        }

        fmpz_set(den, E(pivot_row, pivot_col));
        pivot_row++;
        pivot_col++;
    }

    fmpz_clear(t);

    if (sign < 0)
    {
        fmpz_neg(den, den);
        fmpz_mat_neg(B, B);
    }

    return rank;
}
