/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"

truth_t
_gr_vec_contains(gr_srcptr vec, slong len, gr_srcptr x, gr_ctx_t ctx)
{
    truth_t contains = T_FALSE;

    for (slong i = 0; i < len; i++)
    {
        gr_srcptr y = GR_ENTRY(vec, i, ctx->sizeof_elem);
        contains = truth_or(contains, gr_equal(x, y, ctx));
        if (contains == T_TRUE)
            break;
    }

    return contains;
}
