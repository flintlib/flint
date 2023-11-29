/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

/* Warshall's algorithm */
void
bool_mat_transitive_closure(bool_mat_t dest, const bool_mat_t src)
{
    slong k, i, j, dim;

    if (bool_mat_nrows(dest) != bool_mat_nrows(src) ||
        bool_mat_ncols(dest) != bool_mat_ncols(src))
    {
        flint_throw(FLINT_ERROR, "bool_mat_transitive_closure: incompatible dimensions\n");
    }

    dim = bool_mat_nrows(src);
    if (dim != bool_mat_ncols(src))
    {
        flint_throw(FLINT_ERROR, "bool_mat_transitive_closure: a square matrix is required!\n");
    }

    bool_mat_set(dest, src);

    for (k = 0; k < dim; k++)
    {
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
            {
                if (!bool_mat_get_entry(dest, i, j))
                {
                    bool_mat_set_entry(dest, i, j, (
                                bool_mat_get_entry(dest, i, k) &
                                bool_mat_get_entry(dest, k, j)));
                }
            }
        }
    }
}
