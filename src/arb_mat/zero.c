/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_zero(arb_mat_t mat)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(mat); i++)
        for (j = 0; j < arb_mat_ncols(mat); j++)
            arb_zero(arb_mat_entry(mat, i, j));
}
