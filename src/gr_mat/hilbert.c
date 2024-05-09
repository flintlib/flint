/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"

int
gr_mat_hilbert(gr_mat_t mat, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    /* todo: exploit symmetry */
    for (i = 0; i < R; i++)
    {
        for (j = 0; j < C; j++)
        {
            status |= gr_one(GR_MAT_ENTRY(mat, i, j, sz), ctx);
            status |= gr_div_ui(GR_MAT_ENTRY(mat, i, j, sz),
                GR_MAT_ENTRY(mat, i, j, sz), i + j + 1, ctx);
        }
    }

    return status;
}
