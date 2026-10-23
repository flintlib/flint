/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_randops(gr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx)
{
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);
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
                status |= _gr_vec_add(GR_MAT_ENTRY(mat, j, 0, sz), GR_MAT_ENTRY(mat, j, 0, sz), GR_MAT_ENTRY(mat, i, 0, sz), n, ctx);
            else
                status |= _gr_vec_sub(GR_MAT_ENTRY(mat, j, 0, sz), GR_MAT_ENTRY(mat, j, 0, sz), GR_MAT_ENTRY(mat, i, 0, sz), n, ctx);
        }
        else
        {
            if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < m; k++)
                    status |= add(GR_MAT_ENTRY(mat, k, j, sz), GR_MAT_ENTRY(mat, k, j, sz), GR_MAT_ENTRY(mat, k, i, sz), ctx);
            else
                for (k = 0; k < m; k++)
                    status |= sub(GR_MAT_ENTRY(mat, k, j, sz), GR_MAT_ENTRY(mat, k, j, sz), GR_MAT_ENTRY(mat, k, i, sz), ctx);
        }
    }

    return status;
}
