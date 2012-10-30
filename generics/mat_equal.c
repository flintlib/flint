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

int
elem_mat_equal(const elem_mat_t A, const elem_mat_t B, const ring_t ring)
{
    long i, j;

    if (A->r != B->r || A->c != B->c)
        return 0;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (!elem_equal(elem_mat_entry(A, i, j, ring),
                    elem_mat_entry(B, i, j, ring), RING_PARENT(ring)))
                return 0;
    return 1;
}

