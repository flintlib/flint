/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_hilbert(arb_mat_t mat, slong prec)
{
    slong R, C, i, j;

    R = arb_mat_nrows(mat);
    C = arb_mat_ncols(mat);

    for (i = 0; i < R; i++)
    {
        for (j = 0; j < C; j++)
        {
            arb_one(arb_mat_entry(mat, i, j));
            arb_div_ui(arb_mat_entry(mat, i, j),
                arb_mat_entry(mat, i, j), i + j + 1, prec);
        }
    }
}

