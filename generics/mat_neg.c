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
_elem_mat_neg(elem_ptr * B, const elem_ptr * A, long ar, long ac, const ring_t ring)
{
    long i, j;

    for (i = 0; i < ar; i++)
        for (j = 0; j < ac; j++)
            elem_neg(MAT_INDEX(B, i, j, ring), MAT_SRCINDEX(A, i, j, ring), ring);
}

void
elem_mat_neg(elem_mat_t B, const elem_mat_t A, const ring_t ring)
{
    if (B->r != A->r || B->c != A->c)
    {
        printf("bad dimensions for mat_neg");
        abort();
    }

    _elem_mat_neg(B->rows, A->rows, A->r, A->c, RING_PARENT(ring));
}

