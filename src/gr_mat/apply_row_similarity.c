/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"

int gr_mat_apply_row_similarity(gr_mat_t A, slong r, gr_ptr d, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n = A->r, i, j;
    slong sz = ctx->sizeof_elem;

    if (A->r != A->c || r < 0 || r > A->r)
        return GR_DOMAIN;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < r - 1; j++)
            status |= gr_addmul(GR_MAT_ENTRY(A, i, j, sz), GR_MAT_ENTRY(A, i, r, sz), d, ctx);

        for (j = r + 1; j < n; j++)
            status |= gr_addmul(GR_MAT_ENTRY(A, i, j, sz), GR_MAT_ENTRY(A, i, r, sz), d, ctx);
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < r - 1; j++)
            status |= gr_submul(GR_MAT_ENTRY(A, r, i, sz), GR_MAT_ENTRY(A, j, i, sz), d, ctx);

        for (j = r + 1; j < n; j++)
            status |= gr_submul(GR_MAT_ENTRY(A, r, i, sz), GR_MAT_ENTRY(A, j, i, sz), d, ctx);
    }

    return status;
}
