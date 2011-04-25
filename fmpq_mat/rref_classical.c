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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

static int
_pivot(long * perm, fmpq_mat_t mat, long from_row, long in_column)
{
    long t, j;
    fmpq * u;

    if (!fmpq_is_zero(fmpq_mat_entry(mat, from_row, in_column)))
        return 1;

    for (j = from_row + 1; j < mat->r; j++)
    {
        if (!fmpq_is_zero(fmpq_mat_entry(mat, j, in_column)))
        {
            if (perm)
            {
                t = perm[j];
                perm[j] = perm[from_row];
                perm[from_row] = t;
            }

            u = mat->rows[j];
            mat->rows[j] = mat->rows[from_row];
            mat->rows[from_row] = u; 
            return -1;
        }
    }
    return 0;
}


long
fmpq_mat_rref_classical(long * perm, fmpq_mat_t B, const fmpq_mat_t A)
{
    long m, n, i, j, pivot_row, pivot_col, rank;

    m = A->r;
    n = A->c;

    if (m == 0 || n == 0)
        return 0;

    if (A != B)
        fmpq_mat_set(B, A);

    if (perm)
        for (j = 0; j < m; j++)
            perm[j] = j;

    rank = 0;
    pivot_row = 0;
    pivot_col = 0;

    while (pivot_row < m && pivot_col < n)
    {
        if (!_pivot(perm, B, pivot_row, pivot_col))
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
            if (i == pivot_row)
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

