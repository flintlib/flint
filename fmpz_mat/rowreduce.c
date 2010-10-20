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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

/*
  Helper function for Gaussian elimination. This function searches for k,
  from_row <= k < n, such that a[k][in_column] is nonzero and therefore can
  be used as a pivot element for row reduction. If from_row = k, it
  returns +1. If from_row < k, it swaps row from_row and k and returns -1.
  If no nonzero pivot element is found, it returns 0.
*/
static int
_pivot(fmpz ** rows, long n, long from_row, long in_column)
{
    fmpz * tmp;
    long j;

    if (rows[from_row][in_column] != 0)
        return 1;

    for (j = from_row + 1; j < n; j++)
    {

        if (rows[j][in_column] != 0L)
        {
            tmp = rows[j];
            rows[j] = rows[from_row];
            rows[from_row] = tmp;
            return -1;
        }
    }
    return 0;
}

long _fmpz_mat_rowreduce(fmpz ** a, long m, long n, int options)
{
    long j, k, rank;
    int sign = 1;
    int det_sign;

    long pivot_row;
    long pivot_col;
    long prev_pivot_row;
    long prev_pivot_col;

    fmpz_t d;
    fmpz_init(d);

    rank = 0L;
    prev_pivot_row = -1L;
    prev_pivot_col = -1L;
    pivot_row = 0L;
    pivot_col = 0L;

    while (pivot_row < m && pivot_col < n)
    {
        det_sign = _pivot(a, m, pivot_row, pivot_col);

        if (!det_sign)
        {
            if (options & ROWREDUCE_FAST_ABORT)
            {
                rank = 0L;
                break;
            }
            pivot_col++;
            continue;
        }

        sign *= det_sign;
        rank++;

        if (prev_pivot_row >= 0)
            fmpz_set(d, &a[prev_pivot_row][prev_pivot_col]);

        for (j = (options & ROWREDUCE_FULL ? 0 : pivot_row + 1); j < m; j++)
        {
            if (j == pivot_row)
                continue;
            for (k = pivot_col + 1; k < n; k++)
            {
                fmpz_mul   (&a[j][k], &a[pivot_row][pivot_col], &a[j][k]);
                fmpz_submul(&a[j][k], &a[j][pivot_col], &a[pivot_row][k]);
                if (prev_pivot_row > -1)
                    fmpz_divexact(&a[j][k], &a[j][k], d);
            }

            if (options & ROWREDUCE_CLEAR_LOWER)
                fmpz_zero(&a[j][pivot_col]);
        }

        prev_pivot_row = pivot_row;
        prev_pivot_col = pivot_col;
        pivot_row++;
        pivot_col++;
    }

    fmpz_clear(d);
    return rank * sign;
}
