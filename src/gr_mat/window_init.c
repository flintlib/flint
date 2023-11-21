/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

void
gr_mat_window_init(gr_mat_t window, const gr_mat_t mat,
    slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx)
{
    slong i, sz = ctx->sizeof_elem;
    window->entries = NULL;

    window->rows = flint_malloc((r2 - r1) * sizeof(gr_ptr));

    for (i = 0; i < r2 - r1; i++)
        window->rows[i] = GR_ENTRY(mat->rows[r1 + i], c1, sz);

    window->r = r2 - r1;
    window->c = c2 - c1;
}
