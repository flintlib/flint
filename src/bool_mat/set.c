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
bool_mat_set(bool_mat_t dest, const bool_mat_t src)
{
    slong i, j;

    if (dest == src || bool_mat_is_empty(src))
        return;

    for (i = 0; i < bool_mat_nrows(src); i++)
        for (j = 0; j < bool_mat_ncols(src); j++)
            bool_mat_set_entry(dest, i, j, bool_mat_get_entry(src, i, j));
}
