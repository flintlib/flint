/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"

int
gr_vec_append(gr_vec_t vec, gr_srcptr f, gr_ctx_t ctx)
{
    gr_vec_fit_length(vec, vec->length + 1, ctx);
    vec->length++;
    return gr_set(GR_ENTRY(vec->entries, vec->length - 1, ctx->sizeof_elem), f, ctx);
}
