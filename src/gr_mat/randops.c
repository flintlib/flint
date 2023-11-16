/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_randops(gr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx)
{
    slong c, i, j, k;
    slong m = mat->r;
    slong n = mat->c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (mat->r == 0 || mat->c == 0)
        return GR_SUCCESS;

    for (c = 0; c < count; c++)
    {
        if (n_randint(state, 2))
        {
            if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < n; k++)
                    status |= gr_add(GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, i, k, sz), ctx);
            else
                for (k = 0; k < n; k++)
                    status |= gr_sub(GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, i, k, sz), ctx);
        }
        else
        {
            if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < m; k++)
                    status |= gr_add(GR_MAT_ENTRY(mat, k, j, sz), GR_MAT_ENTRY(mat, k, j, sz), GR_MAT_ENTRY(mat, k, i, sz), ctx);
            else
                for (k = 0; k < m; k++)
                    status |= gr_sub(GR_MAT_ENTRY(mat, k, j, sz), GR_MAT_ENTRY(mat, k, j, sz), GR_MAT_ENTRY(mat, k, i, sz), ctx);
        }
    }

    return status;
}
