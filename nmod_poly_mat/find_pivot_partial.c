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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

len_t
nmod_poly_mat_find_pivot_partial(const nmod_poly_mat_t mat,
                                    len_t start_row, len_t end_row, len_t c)
{
    len_t best_row, best_length, i;

    best_row = start_row;
    best_length = nmod_poly_length(nmod_poly_mat_entry(mat, start_row, c));

    for (i = start_row + 1; i < end_row; i++)
    {
        len_t l;

        l = nmod_poly_length(nmod_poly_mat_entry(mat, i, c));

        if (l != 0 && (best_length == 0 || l <= best_length))
        {
            best_row = i;
            best_length = l;
        }
    }

    if (best_length == 0)
        return -1;

    return best_row;
}
