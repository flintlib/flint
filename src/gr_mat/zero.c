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
gr_mat_zero(gr_mat_t res, gr_ctx_t ctx)
{
    int status;
    slong i, r, c;
    slong sz = ctx->sizeof_elem;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    status = GR_SUCCESS;
    for (i = 0; i < r; i++)
    {
        status |= _gr_vec_zero(GR_MAT_ENTRY(res, i, 0, sz), c, ctx);
    }

    return status;
}
