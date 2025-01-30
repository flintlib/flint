/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "long_extras.h"
#include "gr.h"
#include "gr_mat.h"

void
gr_mat_init(gr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)
{
    mat->entries = NULL;
    mat->r = rows;
    mat->c = cols;
    mat->stride = cols;

    if (rows != 0 && cols != 0)
    {
        slong sz = ctx->sizeof_elem;
        slong num;
        int of;

        of = z_mul_checked(&num, rows, cols);

        if (of)
            flint_throw(FLINT_ERROR, "Overflow creating a %wd x %wd object\n", rows, cols);

        mat->entries = flint_malloc(num * sz);
        _gr_vec_init(mat->entries, rows * cols, ctx);
    }
}
