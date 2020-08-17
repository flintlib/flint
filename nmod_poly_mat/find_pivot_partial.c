/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

slong
nmod_poly_mat_find_pivot_partial(const nmod_poly_mat_t mat,
                                    slong start_row, slong end_row, slong c)
{
    slong best_row, best_length, i;

    best_row = start_row;
    best_length = nmod_poly_length(nmod_poly_mat_entry(mat, start_row, c));

    for (i = start_row + 1; i < end_row; i++)
    {
        slong l;

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
