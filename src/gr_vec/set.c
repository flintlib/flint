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
gr_vec_set(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (res != src)
    {
        gr_vec_set_length(res, src->length, ctx);
        status = _gr_vec_set(res->entries, src->entries, src->length, ctx);
    }

    return status;
}
