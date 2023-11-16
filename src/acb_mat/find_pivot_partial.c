/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_cmpabs_approx(const acb_t x, const acb_t y)
{
    const arf_struct *xm, *ym;

    if (arf_cmpabs(arb_midref(acb_realref(x)), arb_midref(acb_imagref(x))) >= 0)
        xm = arb_midref(acb_realref(x));
    else
        xm = arb_midref(acb_imagref(x));

    if (arf_cmpabs(arb_midref(acb_realref(y)), arb_midref(acb_imagref(y))) >= 0)
        ym = arb_midref(acb_realref(y));
    else
        ym = arb_midref(acb_imagref(y));

    return arf_cmpabs(xm, ym);
}

slong
acb_mat_find_pivot_partial(const acb_mat_t mat,
                                    slong start_row, slong end_row, slong c)
{
    slong best_row, i;

    best_row = -1;

    for (i = start_row; i < end_row; i++)
    {
        if (!acb_contains_zero(acb_mat_entry(mat, i, c)))
        {
            if (best_row == -1)
            {
                best_row = i;
            }
            /* todo: should take the radius into account */
            else if (acb_cmpabs_approx(acb_mat_entry(mat, i, c), acb_mat_entry(mat, best_row, c)) > 0)
            {
                best_row = i;
            }
        }
    }
    return best_row;
}

