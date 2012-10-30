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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

long
elem_mat_rref(elem_mat_t R, elem_ptr den, const elem_mat_t A, const ring_t ring)
{
    long i, j, k, m, n, rank;
    long *pivots, *nonpivots;
    const ring_struct * ering = RING_PARENT(ring);

    rank = elem_mat_fflu(R, den, NULL, A, 0, ring);
    m = elem_mat_nrows(R, ring);
    n = elem_mat_ncols(R, ring);

    /* clear bottom */
    for (i = rank; i < m; i++)
        for (j = 0; j < n; j++)
            elem_zero(elem_mat_entry(R, i, j, ring), ering);

    /* Convert row echelon form to reduced row echelon form */
    if (rank > 1)
    {
        elem_ptr tmp, tmp2;

        ELEM_TMP_INIT(tmp, ering);
        ELEM_TMP_INIT(tmp2, ering);

        pivots = flint_malloc(sizeof(long) * n);
        nonpivots = pivots + rank;

        /* find pivot positions */
        for (i = j = k = 0; i < rank; i++)
        {
            while (elem_is_zero(elem_mat_entry(R, i, j, ring), ering))
            {
                nonpivots[k] = j;
                k++;
                j++;
            }
            pivots[i] = j;
            j++;
        }
        while (k < n - rank)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        for (k = 0; k < n - rank; k++)
        {
            for (i = rank - 2; i >= 0; i--)
            {
                elem_mul(tmp, den, elem_mat_entry(R, i, nonpivots[k], ring), ering);

                for (j = i + 1; j < rank; j++)
                {
                    elem_mul(tmp2, elem_mat_entry(R, i, pivots[j], ring), elem_mat_entry(R, j, nonpivots[k], ring), ering);
                    elem_sub(tmp, tmp, tmp2, ering);
                }

                elem_divexact(elem_mat_entry(R, i, nonpivots[k], ring), tmp, elem_mat_entry(R, i, pivots[i], ring), ering);
            }
        }

        /* clear pivot columns */
        for (i = 0; i < rank; i++)
        {
            for (j = 0; j < rank; j++)
            {
                if (i == j)
                    elem_set(elem_mat_entry(R, j, pivots[i], ring), den, ering);
                else
                    elem_zero(elem_mat_entry(R, j, pivots[i], ring), ering);
            }
        }

        flint_free(pivots);
        ELEM_TMP_CLEAR(tmp, ering);
        ELEM_TMP_CLEAR(tmp2, ering);
    }

    return rank;
}


