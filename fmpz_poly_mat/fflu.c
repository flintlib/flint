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

static __inline__ void
fmpz_poly_mat_swap_rows(fmpz_poly_mat_t mat, len_t * perm, len_t r, len_t s)
{
    if (r != s)
    {
        fmpz_poly_struct * u;
        len_t t;

        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u; 
    }
}

len_t
fmpz_poly_mat_fflu(fmpz_poly_mat_t B, fmpz_poly_t den, len_t * perm,
    const fmpz_poly_mat_t A, int rank_check)
{
    fmpz_poly_t t;
    len_t m, n, j, k, rank, r, pivot_row, pivot_col;

    if (fmpz_poly_mat_is_empty(A))
    {
        fmpz_poly_one(den);
        return 0;
    }

    fmpz_poly_mat_set(B, A);
    m = B->r;
    n = B->c;
    rank = pivot_row = pivot_col = 0;

    fmpz_poly_init(t);

    while (pivot_row < m && pivot_col < n)
    {
        r = fmpz_poly_mat_find_pivot_partial(B, pivot_row, m, pivot_col);

        if (r == -1)
        {
            if (rank_check)
            {
                fmpz_poly_zero(den);
                rank = 0;
                break;
            }
            pivot_col++;
            continue;
        }
        else if (r != pivot_row)
            fmpz_poly_mat_swap_rows(B, perm, pivot_row, r);

        rank++;

        for (j = pivot_row + 1; j < m; j++)
        {
            for (k = pivot_col + 1; k < n; k++)
            {
                fmpz_poly_mul(E(j, k), E(j, k), E(pivot_row, pivot_col));
                fmpz_poly_mul(t, E(j, pivot_col), E(pivot_row, k));
                fmpz_poly_sub(E(j, k), E(j, k), t);
                if (pivot_row > 0)
                    fmpz_poly_div(E(j, k), E(j, k), den);
            }
        }

        fmpz_poly_set(den, E(pivot_row, pivot_col));
        pivot_row++;
        pivot_col++;
    }

    fmpz_poly_clear(t);
    return rank;
}
