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

int
gr_mat_set_gr_mat_other(gr_mat_t res, const gr_mat_t mat, gr_ctx_t mat_ctx, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    slong mat_sz = mat_ctx->sizeof_elem;
    int status = GR_SUCCESS;

    R = gr_mat_nrows(mat, mat_ctx);
    C = gr_mat_ncols(mat, mat_ctx);

    if (R != gr_mat_nrows(res, ctx) || C != gr_mat_ncols(res, ctx))
        return GR_DOMAIN;

    for (i = 0; i < R; i++)
        for (j = 0; j < C && status == GR_SUCCESS; j++)
            status |= gr_set_other(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, mat_sz), mat_ctx, ctx);

    return status;
}