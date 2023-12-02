/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

void
bool_mat_directed_path(bool_mat_t mat)
{
    slong i;

    if (!bool_mat_is_square(mat))
    {
        flint_throw(FLINT_ERROR, "bool_mat_directed_path: a square matrix is required!\n");
    }

    if (bool_mat_is_empty(mat))
        return;

    bool_mat_zero(mat);
    for (i = 0; i < bool_mat_nrows(mat) - 1; i++)
    {
        bool_mat_set_entry(mat, i, i+1, 1);
    }
}
