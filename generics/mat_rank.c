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

long
elem_mat_rank(const elem_mat_t mat, const ring_t ring)
{
    elem_mat_t tmp;
    elem_ptr den;
    long m, n, rank;

    m = elem_mat_nrows(mat, ring);
    n = elem_mat_ncols(mat, ring);

    if (m == 0 || n == 0)
        return 0;

    elem_mat_init(tmp, m, n, ring);
    ELEM_TMP_INIT(den, RING_PARENT(ring));
    rank = elem_mat_fflu(tmp, den, NULL, mat, 0, ring);
    ELEM_TMP_CLEAR(den, RING_PARENT(ring));
    elem_mat_clear(tmp, ring);

    return rank;
}

