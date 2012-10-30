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
elem_mat_scalar_mul(elem_mat_t res, const elem_mat_t mat, elem_srcptr c, const ring_t ring)
{
    long i, j;

    for (i = 0; i < elem_mat_nrows(mat, ring); i++)
        for (j = 0; j < elem_mat_ncols(mat, ring); j++)
            elem_mul(elem_mat_entry(res, i, j, ring), elem_mat_entry(mat, i, j, ring), c, RING_PARENT(ring));
}

