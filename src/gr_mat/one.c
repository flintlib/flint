/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_one(gr_mat_t res, gr_ctx_t ctx)
{
    int status;
    slong i, r, c, sz;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);
    sz = ctx->sizeof_elem;

    status = gr_mat_zero(res, ctx);

    if (r > 0 && c > 0)
    {
        for (i = 0; i < FLINT_MIN(r, c); i++)
            status |= gr_one(GR_MAT_ENTRY(res, i, i, sz), ctx);
    }

    return status;
}
