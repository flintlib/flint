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
gr_vec_fit_length(gr_vec_t vec, slong len, gr_ctx_t ctx)
{
    slong alloc = vec->alloc;

    if (len > alloc)
    {
        slong sz = ctx->sizeof_elem;

        if (len < 2 * alloc)
            len = 2 * alloc;

        vec->entries = flint_realloc(vec->entries, len * sz);
        _gr_vec_init(GR_ENTRY(vec->entries, vec->alloc, sz), len - alloc, ctx);
        vec->alloc = len;
    }
}
