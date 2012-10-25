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
_elem_mat_transpose_square_inplace(elem_ptr * A, long ar, long ac, const ring_t ring)
{
    long i, j;

    for (i = 0; i < ar - 1; i++)
        for (j = i + 1; j < ac; j++)
            elem_swap(MAT_INDEX(A, i, j, ring), MAT_INDEX(A, j, i, ring), ring);
}

void
_elem_mat_transpose(elem_ptr * B, elem_ptr * A, long ar, long ac, const ring_t ring)
{
    long i, j;

    for (i = 0; i < ac; i++)
        for (j = 0; j < ar; j++)
            elem_set(MAT_INDEX(B, i, j, ring), MAT_INDEX(A, j, i, ring), ring);
}

void
elem_mat_transpose(elem_mat_t B, const elem_mat_t A, const ring_t ring)
{
    if (B->r != A->c || B->c != A->r)
    {
        printf("bad dimensions for mat_transpose\n");
        abort();
    }

    if (A == B)
        _elem_mat_transpose_square_inplace(B->rows, B->r, B->c, RING_PARENT(ring));
    else
        _elem_mat_transpose(B->rows, A->rows, A->r, A->c, RING_PARENT(ring));
}

