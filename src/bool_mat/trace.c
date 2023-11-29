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
bool_mat_trace(const bool_mat_t mat)
{
    slong i;

    if (!bool_mat_is_square(mat))
    {
        flint_throw(FLINT_ERROR, "bool_mat_trace: a square matrix is required!\n");
    }

    if (bool_mat_is_empty(mat))
        return 0;

    for (i = 0; i < bool_mat_nrows(mat); i++)
        if (bool_mat_get_entry(mat, i, i))
            return 1;

    return 0;
}
