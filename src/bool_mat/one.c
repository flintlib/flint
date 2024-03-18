/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

void
bool_mat_one(bool_mat_t mat)
{
    slong i, j;

    for (i = 0; i < bool_mat_nrows(mat); i++)
        for (j = 0; j < bool_mat_ncols(mat); j++)
            bool_mat_set_entry(mat, i, j, i == j);
}
