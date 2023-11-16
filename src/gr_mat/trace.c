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
gr_mat_trace(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong i, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = mat->r;
    c = mat->c;

    if (r != c)
        return GR_DOMAIN;

    if (r == 0)
        return gr_zero(res, ctx);

    if (r == 1)
        return gr_set(res, GR_MAT_ENTRY(mat, 0, 0, sz), ctx);

    status |= gr_add(res, GR_MAT_ENTRY(mat, 0, 0, sz),
        GR_MAT_ENTRY(mat, 1, 1, sz), ctx);

    for (i = 2; i < r; i++)
        status |= gr_add(res, res, GR_MAT_ENTRY(mat, i, i, sz), ctx);

    return status;
}
