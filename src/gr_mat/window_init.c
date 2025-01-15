/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

void
gr_mat_window_init(gr_mat_t window, const gr_mat_t mat,
    slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;

    FLINT_ASSERT(r1 >= 0 && r1 <= r2 && r2 <= mat->r);
    FLINT_ASSERT(c2 >= 0 && c1 <= c2 && c2 <= mat->c);

    window->entries = GR_MAT_ENTRY(mat, r1, c1, sz);
    window->r = r2 - r1;
    window->c = c2 - c1;
    window->stride = mat->stride;
}
