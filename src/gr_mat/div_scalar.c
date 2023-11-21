/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/* todo: use a vector function; preinvert when appropriate */
int
gr_mat_div_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
{
    slong i, j, r, c;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (c != 0)
        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
                status |= gr_div(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), x, ctx);

    return status;
}
