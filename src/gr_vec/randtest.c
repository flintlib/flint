/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"

int
_gr_vec_randtest(gr_ptr res, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;

    status = GR_SUCCESS;
    for (i = 0; i < len; i++)
    {
        if (n_randint(state, 2))
            status |= gr_zero(GR_ENTRY(res, i, sz), ctx);
        else
            status |= gr_randtest(GR_ENTRY(res, i, sz), state, ctx);
    }

    return status;
}
