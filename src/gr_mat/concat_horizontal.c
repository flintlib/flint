/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_concat_horizontal(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong c2 = mat2->c;
    slong sz = ctx->sizeof_elem;

    if (mat1->r != mat2->r || res->c != mat1->c + mat2->c)
        return GR_DOMAIN;

    for (i = 0; i < r1; i++)
    {
        if (c1 > 0)
            status |= _gr_vec_set(res->rows[i], mat1->rows[i], c1, ctx);
        if (c2 > 0)
            status |= _gr_vec_set(GR_ENTRY(res->rows[i], c1, sz), mat2->rows[i], c2, ctx);
    }

    return status;
}
