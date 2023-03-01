/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
arb_mat_ne(const arb_mat_t mat1, const arb_mat_t mat2)
{
    slong i, j;

    if ((arb_mat_nrows(mat1) != arb_mat_nrows(mat2)) ||
        (arb_mat_ncols(mat1) != arb_mat_ncols(mat2)))
        return 1;

    for (i = 0; i < arb_mat_nrows(mat1); i++)
        for (j = 0; j < arb_mat_ncols(mat1); j++)
            if (arb_ne(arb_mat_entry(mat1, i, j),
                arb_mat_entry(mat2, i, j)))
                return 1;

    return 0;
}
