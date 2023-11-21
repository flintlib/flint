/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

slong
arb_mat_find_pivot_partial(const arb_mat_t mat,
                                    slong start_row, slong end_row, slong c)
{
    slong best_row, i;

    best_row = -1;

    for (i = start_row; i < end_row; i++)
    {
        if (!arb_contains_zero(arb_mat_entry(mat, i, c)))
        {
            if (best_row == -1)
            {
                best_row = i;
            }
            /* todo: should take the radius into account */
            else if (arf_cmpabs(arb_midref(arb_mat_entry(mat, i, c)),
                                arb_midref(arb_mat_entry(mat, best_row, c))) > 0)
            {
                best_row = i;
            }
        }
    }

    return best_row;
}
