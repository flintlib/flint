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
gr_mat_init(gr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)
{
    slong i, sz;

    sz = ctx->sizeof_elem;

    if (rows != 0)
        mat->rows = flint_malloc(rows * sizeof(gr_ptr));
    else
        mat->rows = NULL;

    if (rows != 0 && cols != 0)
    {
        mat->entries = (gr_ptr) flint_malloc(flint_mul_sizes(rows, cols) * sz);

        _gr_vec_init(mat->entries, rows * cols, ctx);

        for (i = 0; i < rows; i++)
            mat->rows[i] = GR_ENTRY(mat->entries, i * cols, sz);
    }
    else
    {
        mat->entries = NULL;
        for (i = 0; i < rows; i++)
            mat->rows[i] = NULL;
    }

    mat->r = rows;
    mat->c = cols;
}
