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
bool_mat_mul_entrywise(bool_mat_t C, const bool_mat_t A, const bool_mat_t B)
{
    slong i, j;

    if (bool_mat_nrows(A) != bool_mat_nrows(B) ||
        bool_mat_ncols(A) != bool_mat_ncols(B))
    {
        flint_throw(FLINT_ERROR, "bool_mat_mul_entrywise: incompatible dimensions\n");
    }

    for (i = 0; i < bool_mat_nrows(A); i++)
    {
        for (j = 0; j < bool_mat_ncols(A); j++)
        {
            bool_mat_set_entry(C, i, j, (bool_mat_get_entry(A, i, j) &
                                         bool_mat_get_entry(B, i, j)));
        }
    }
}
