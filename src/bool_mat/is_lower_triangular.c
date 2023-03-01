/*
    Copyright (C) 2016 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

int
bool_mat_is_lower_triangular(const bool_mat_t mat)
{
    slong i, j;

    if (bool_mat_is_empty(mat))
        return 1;

    for (j = 0; j < bool_mat_ncols(mat); j++)
        for (i = 0; i < bool_mat_nrows(mat) && i < j; i++)
            if (bool_mat_get_entry(mat, i, j))
                return 0;

    return 1;
}
