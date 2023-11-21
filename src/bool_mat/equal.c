/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

int
bool_mat_equal(const bool_mat_t mat1, const bool_mat_t mat2)
{
    slong i, j;

    if ((bool_mat_nrows(mat1) != bool_mat_nrows(mat2)) ||
        (bool_mat_ncols(mat1) != bool_mat_ncols(mat2)))
        return 0;

    for (i = 0; i < bool_mat_nrows(mat1); i++)
        for (j = 0; j < bool_mat_ncols(mat1); j++)
            if (bool_mat_get_entry(mat1, i, j) !=
                bool_mat_get_entry(mat2, i, j))
            {
                return 0;
            }

    return 1;
}
