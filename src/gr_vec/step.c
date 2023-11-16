/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"

/* todo: with floats, prefer multiplying to adding */
int _gr_vec_step(gr_ptr vec, gr_srcptr start, gr_srcptr step, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status = GR_SUCCESS;
    slong i, sz = ctx->sizeof_elem;;

    if (len <= 0)
        return GR_SUCCESS;

    status |= gr_set(vec, start, ctx);
    for (i = 1; i < len; i++)
        status |= add(GR_ENTRY(vec, i, sz), GR_ENTRY(vec, i - 1, sz), step, ctx);

    return status;
}
