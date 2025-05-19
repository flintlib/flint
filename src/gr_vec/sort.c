/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"

typedef struct
{
    gr_ctx_struct * gr_ctx;
    int status;
}
cmp_ctx_struct;

static int
cmp(const void * x, const void * y, void * ctx)
{
    int res;
    cmp_ctx_struct * cmp_ctx = ctx;
    int status = gr_cmp(&res, x, y, cmp_ctx->gr_ctx);
    if (status != GR_SUCCESS)
        res = 0;
    cmp_ctx->status |= status;
    return res;
}

int
_gr_vec_sort(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    cmp_ctx_struct cmp_ctx = { .gr_ctx = ctx, .status = GR_SUCCESS };
    flint_sort(vec, len, ctx->sizeof_elem, cmp, &cmp_ctx);
    return cmp_ctx.status;
}

int
gr_vec_sort(gr_vec_t dest, const gr_vec_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (dest != src)
        status |= gr_vec_set(dest, src, ctx);

    status |= _gr_vec_sort(dest->entries, dest->length, ctx);

    return status;
}
