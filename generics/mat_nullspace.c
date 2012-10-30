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
elem_mat_nullspace(elem_mat_t res, const elem_mat_t mat, const ring_t ring)
{
    long i, j, k, m, n, rank, nullity;
    long * pivots;
    long * nonpivots;
    const ring_struct * ering = RING_PARENT(ring); 

    elem_mat_t tmp;
    elem_ptr den;

    m = elem_mat_nrows(mat, ring);
    n = elem_mat_ncols(mat, ring);

    ELEM_TMP_INIT(den, ering);
    elem_mat_init(tmp, m, n, ring);
    elem_mat_set(tmp, mat, ring);

    rank = elem_mat_rref(tmp, den, tmp, ring);
    nullity = n - rank;

    elem_mat_zero(res, ring);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            elem_one(elem_mat_entry(res, i, i, ring), ering);
    }
    else if (nullity)
    {
        pivots = flint_malloc(rank * sizeof(long));
        nonpivots = flint_malloc(nullity * sizeof(long));

        for (i = j = k = 0; i < rank; i++)
        {
            while (elem_is_zero(elem_mat_entry(tmp, i, j, ring), ering))
            {
                nonpivots[k] = j;
                k++;
                j++;
            }
            pivots[i] = j;
            j++;
        }
        while (k < nullity)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
                elem_set(elem_mat_entry(res, pivots[j], i, ring),
                    elem_mat_entry(tmp, j, nonpivots[i], ring), ering);

            elem_neg(elem_mat_entry(res, nonpivots[i], i, ring), den, ering);
        }

        flint_free(pivots);
        flint_free(nonpivots);
    }

    ELEM_TMP_CLEAR(den, ering);
    elem_mat_clear(tmp, ring);
    return nullity;
}

