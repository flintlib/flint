/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_ones(acb_mat_t mat)
{
    slong R, C, i, j;

    R = acb_mat_nrows(mat);
    C = acb_mat_ncols(mat);

    for (i = 0; i < R; i++)
        for (j = 0; j < C; j++)
            acb_one(acb_mat_entry(mat, i, j));
}

