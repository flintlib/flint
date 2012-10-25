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
elem_mat_init(elem_mat_t mat, long rows, long cols, const ring_t ring)
{
    if (rows && cols)
    {
        long i, size = RING_PARENT(ring)->size;

        mat->entries = _elem_vec_init(rows * cols, RING_PARENT(ring));
        mat->rows = flint_malloc(rows * sizeof(elem_ptr));

        for (i = 0; i < rows; i++)
            mat->rows[i] = INDEX(mat->entries, i * cols, size);
    }
    else
        mat->entries = NULL;

    mat->r = rows;
    mat->c = cols;
}

