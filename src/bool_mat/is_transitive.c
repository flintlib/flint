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
bool_mat_is_transitive(const bool_mat_t mat)
{
    slong n, i, j, k;

    if (!bool_mat_is_square(mat))
    {
        flint_throw(FLINT_ERROR, "bool_mat_is_transitive: a square matrix is required!\n");
    }

    if (bool_mat_is_empty(mat))
        return 1;

    n = bool_mat_nrows(mat);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                if (bool_mat_get_entry(mat, i, j) &&
                    bool_mat_get_entry(mat, j, k) &&
                    !bool_mat_get_entry(mat, i, k))
                {
                    return 0;
                }

    return 1;
}
