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
_elem_mat_print(const elem_ptr * rows, long r, long c, const ring_t ring)
{
    if (r == 0 && c == 0)
    {
        printf("{{}}");
    }
    else
    {
        long i, j;

        printf("{");

        for (i = 0; i < r; i++)
        {
            printf("{");

            for (j = 0; j < c; j++)
            {
                elem_print(MAT_SRCINDEX(rows, i, j, ring), ring);

                if (j < c - 1)
                    printf(", ");
            }

            printf("}");

            if (i < r - 1)
                printf(",\n");
        }

        printf("}");

    }
}

void
elem_mat_print(const elem_mat_t mat, const ring_t ring)
{
    _elem_mat_print(mat->rows, mat->r, mat->c, RING_PARENT(ring));
}

