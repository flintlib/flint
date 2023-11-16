/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_randops(ca_mat_t mat, flint_rand_t state, slong count, ca_ctx_t ctx)
{
    slong c, i, j, k;
    slong m = mat->r;
    slong n = mat->c;

    if (mat->r == 0 || mat->c == 0)
        return;

    for (c = 0; c < count; c++)
    {
        if (n_randint(state, 2))
        {
            if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < n; k++)
                    ca_add(ca_mat_entry(mat, j, k), ca_mat_entry(mat, j, k), ca_mat_entry(mat, i, k), ctx);
            else
                for (k = 0; k < n; k++)
                    ca_sub(ca_mat_entry(mat, j, k), ca_mat_entry(mat, j, k), ca_mat_entry(mat, i, k), ctx);
        }
        else
        {
            if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < m; k++)
                    ca_add(ca_mat_entry(mat, k, j), ca_mat_entry(mat, k, j), ca_mat_entry(mat, k, i), ctx);
            else
                for (k = 0; k < m; k++)
                    ca_sub(ca_mat_entry(mat, k, j), ca_mat_entry(mat, k, j), ca_mat_entry(mat, k, i), ctx);
        }
    }
}
