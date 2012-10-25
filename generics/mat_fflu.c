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

static __inline__ void
_elem_mat_swap_rows(elem_ptr * mat, long r, long s, const ring_t ring)
{
    elem_ptr u = mat[s];
    mat[s] = mat[r];
    mat[r] = u;
}

static __inline__ void
_perm_swap_rows(long * perm, long r, long s)
{
    long u = perm[s];
    perm[s] = perm[r];
    perm[r] = u;
}


static __inline__ long
_elem_mat_find_pivot_any(const elem_ptr * mat, long start_row, long end_row, long c, const ring_t ring)
{
    long r;

    for (r = start_row; r < end_row; r++)
    {
        if (!elem_is_zero(MAT_SRCINDEX(mat, r, c, ring), ring))
            return r;
    }

    return -1;
}


long
_elem_mat_fflu(elem_ptr * mat, elem_ptr den, long * perm, long m, long n, int rank_check, const ring_t ring)
{
    long j, k, rank, r, pivot_row, pivot_col;
    elem_ptr t;

    ELEM_TMP_INIT(t, ring);
    rank = pivot_row = pivot_col = 0;

    while (pivot_row < m && pivot_col < n)
    {
        r = _elem_mat_find_pivot_any(mat, pivot_row, m, pivot_col, ring);

        if (r == -1)
        {
            if (rank_check)
            {
                elem_zero(den, ring);
                rank = 0;
                break;
            }
            pivot_col++;
            continue;
        }
        else if (r != pivot_row)
        {
            _elem_mat_swap_rows(mat, pivot_row, r, ring);
            if (perm != NULL)
                _perm_swap_rows(perm, pivot_row, r);
        }

        rank++;

        for (j = pivot_row + 1; j < m; j++)
        {
            for (k = pivot_col + 1; k < n; k++)
            {
                elem_mul(MAT_INDEX(mat, j, k, ring), MAT_SRCINDEX(mat, j, k, ring), MAT_SRCINDEX(mat, pivot_row, pivot_col, ring), ring);
                elem_mul(t, MAT_SRCINDEX(mat, j, pivot_col, ring), MAT_SRCINDEX(mat, pivot_row, k, ring), ring);
                if (pivot_row > 0)
                {
                    elem_sub(t, MAT_INDEX(mat, j, k, ring), t, ring);
                    elem_divexact(MAT_INDEX(mat, j, k, ring), t, den, ring);
                }
                else
                {
                    elem_sub(MAT_INDEX(mat, j, k, ring), MAT_INDEX(mat, j, k, ring), t, ring);
                }
            }
        }

        elem_set(den, MAT_SRCINDEX(mat, pivot_row, pivot_col, ring), ring);
        pivot_row++;
        pivot_col++;
    }

    ELEM_TMP_CLEAR(t, ring);
    return rank;
}

long
elem_mat_fflu(elem_mat_t B, elem_ptr den, long * perm, const elem_mat_t A, int rank_check, const ring_t ring)
{
    long rank;

    if (elem_mat_is_empty(A, ring))
    {
        elem_one(den, RING_PARENT(ring));
        return 0;
    }

    if (A != B)
        elem_mat_set(B, A, ring);

    rank = _elem_mat_fflu(B->rows, den, perm, B->r, B->c, rank_check, RING_PARENT(ring));
    return rank;
}

