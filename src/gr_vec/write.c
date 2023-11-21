/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gr_vec.h"

int
_gr_vec_write(gr_stream_t out, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, sz = ctx->sizeof_elem;

    gr_stream_write(out, "[");
    for (i = 0; i < len; i++)
    {
        status |= gr_write(out, GR_ENTRY(vec, i, sz), ctx);
        if (i < len - 1)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "]");
    return status;
}

int
gr_vec_write(gr_stream_t out, const gr_vec_t vec, gr_ctx_t ctx)
{
    return _gr_vec_write(out, vec->entries, vec->length, ctx);
}

int
_gr_vec_print(gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return _gr_vec_write(out, vec, len, ctx);
}


int
gr_vec_print(const gr_vec_t vec, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_vec_write(out, vec, ctx);
}
