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
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

#define E(j,k) fmpz_poly_mat_entry(B,j,k)

long
fmpz_poly_mat_rowreduce(long * perm, fmpz_poly_mat_t B, fmpz_poly_t den,
    const fmpz_poly_mat_t A, int options)
{
    fmpz_poly_t t;
    long m, n, j, k, rank, pivot_row, pivot_col;
    int sign, temp_sign;

    if (fmpz_poly_mat_is_empty(A))
    {
        fmpz_poly_set_ui(den, 1UL);
        return 0;
    }

    fmpz_poly_mat_set(B, A);
    m = B->r;
    n = B->c;
    rank = pivot_row = pivot_col = 0;
    sign = 1;

    fmpz_poly_init(t);

    while (pivot_row < m && pivot_col < n)
    {
        temp_sign = fmpz_poly_mat_pivot(perm, B, pivot_row, pivot_col);

        if (temp_sign == 0)
        {
            if (options & ROWREDUCE_FAST_ABORT)
            {
                rank = 0;
                break;
            }
            pivot_col++;
            continue;
        }

        sign *= temp_sign;
        rank++;

        for (j = pivot_row + 1; j < m; j++)
        {
            if (j == pivot_row)
                continue;

            for (k = (!(options & ROWREDUCE_FULL) || j > pivot_row) ?
                pivot_col + 1 : j; k < n; k++)
            {
                if (k == pivot_col)
                    continue;

                fmpz_poly_mul(E(j, k), E(j, k), E(pivot_row, pivot_col));
                fmpz_poly_mul(t, E(j, pivot_col), E(pivot_row, k));
                fmpz_poly_sub(E(j, k), E(j, k), t);

                if (pivot_row > 0)
                    fmpz_poly_div(E(j, k), E(j, k), den);
            }

            if (options & ROWREDUCE_CLEAR_LOWER || options & ROWREDUCE_FULL)
                fmpz_poly_zero(E(j, pivot_col));
        }

        fmpz_poly_set(den, E(pivot_row, pivot_col));
        pivot_row++;
        pivot_col++;
    }

    fmpz_poly_clear(t);

    if (sign < 0)
    {
        fmpz_poly_neg(den, den);
        fmpz_poly_mat_neg(B, B);
    }

    return rank;
}
