/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"

void
gr_vec_set_length(gr_vec_t vec, slong len, gr_ctx_t ctx)
{
    if (vec->length > len)
    {
        /* potentially free no longer used elements */
        GR_MUST_SUCCEED(_gr_vec_zero(GR_ENTRY(vec->entries, len, ctx->sizeof_elem), vec->length - len, ctx));
    }
    else if (vec->length < len)
    {
        gr_vec_fit_length(vec, len, ctx);
        GR_MUST_SUCCEED(_gr_vec_zero(GR_ENTRY(vec->entries, vec->length, ctx->sizeof_elem), len - vec->length, ctx));
    }

    vec->length = len;
}
