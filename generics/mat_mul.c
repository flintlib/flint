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

void
_elem_mat_mul(elem_ptr * C, const elem_ptr * A, long ar, long ac,
                const elem_ptr * B, long bc, const ring_t ring)
{
    long i, j, k;
    elem_ptr t;

    ELEM_TMP_INIT(t, ring);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            elem_mul(MAT_INDEX(C, i, j, ring),
                MAT_SRCINDEX(A, i, 0, ring), MAT_SRCINDEX(B, 0, j, ring), ring);

            for (k = 1; k < ac; k++)
            {
                elem_mul(t, MAT_SRCINDEX(A, i, k, ring), MAT_SRCINDEX(B, k, j, ring), ring);
                elem_add(MAT_INDEX(C, i, j, ring), MAT_SRCINDEX(C, i, j, ring), t, ring);
            }
        }
    }

    ELEM_TMP_CLEAR(t, ring);
}

void
elem_mat_mul(elem_mat_t C, const elem_mat_t A, const elem_mat_t B, const ring_t ring)
{
    if (A->c != B->r || C->r != A->r || C->c != B->c)
    {
        printf("bad dimensions for mat_mul");
        abort();
    }

    if (A->c == 0)
    {
        elem_mat_zero(C, ring);
        return;
    }

    if (A == C || B == C)
    {
        elem_mat_t T;
        elem_mat_init(T, C->r, C->c, ring);
        _elem_mat_mul(T->rows, A->rows, A->r, A->c, B->rows, B->c, RING_PARENT(ring));
        elem_mat_swap(C, T, ring);
        elem_mat_clear(T, ring);
    }
    else
    {
        _elem_mat_mul(C->rows, A->rows, A->r, A->c, B->rows, B->c, RING_PARENT(ring));
    }
}

