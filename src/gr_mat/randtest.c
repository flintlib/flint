/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

/* todo: want to try different algorithms here */
int
gr_mat_randtest(gr_mat_t mat, flint_rand_t state, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;

    r = gr_mat_nrows(mat, ctx);
    c = gr_mat_ncols(mat, ctx);

    if (n_randint(state, 10) == 0)
    {
        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
                status |= gr_randtest(GR_MAT_ENTRY(mat, i, j, sz), state, ctx);
    }
    else
    {
        for (i = 0; i < r; i++)
            status |= _gr_vec_randtest(mat->rows[i], state, c, ctx);
    }

    return status;
}
