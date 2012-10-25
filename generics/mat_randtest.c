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
_elem_mat_randtest(elem_ptr * mat,
    flint_rand_t state, long r, long c, const long * size, const ring_t ring)
{
    long i, j;

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            elem_randtest(MAT_INDEX(mat, i, j, ring), state, size, ring);
}

void
elem_mat_randtest(elem_mat_t mat, flint_rand_t state, const long * size, const ring_t ring)
{
    _elem_mat_randtest(mat->rows, state, mat->r, mat->c, size, RING_PARENT(ring));
}

