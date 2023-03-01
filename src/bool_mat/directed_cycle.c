/*
    Copyright (C) 2016 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

void
bool_mat_directed_cycle(bool_mat_t mat)
{
    if (!bool_mat_is_square(mat))
    {
        flint_printf("bool_mat_directed_cycle: a square matrix is required!\n");
        flint_abort();
    }

    if (bool_mat_is_empty(mat))
        return;

    bool_mat_directed_path(mat);
    bool_mat_set_entry(mat, bool_mat_nrows(mat) - 1, 0, 1);
}
