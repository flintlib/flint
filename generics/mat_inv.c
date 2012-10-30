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
#include "perm.h"

int
elem_mat_inv(elem_mat_t X, elem_ptr den, const elem_mat_t A, const ring_t ring)
{
    if (X == A)
    {
        int r;
        elem_mat_t T;
        elem_mat_init(T, elem_mat_nrows(A, ring), elem_mat_ncols(A, ring), ring);
        r = elem_mat_inv(T, den, A, ring);
        elem_mat_swap(T, X, ring);
        elem_mat_clear(T, ring);
        return r;
    }

    elem_mat_one(X, ring);
    return elem_mat_solve(X, den, A, X, ring);
}

