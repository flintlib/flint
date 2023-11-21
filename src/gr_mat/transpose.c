/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_transpose(gr_mat_t B, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (B->r != A->c || B->c != A->r)
        return GR_DOMAIN;

    if (A->r == 0 || A->c == 0)
        return GR_SUCCESS;

    if (A == B)  /* In-place, guaranteed to be square */
    {
        for (i = 0; i < A->r - 1; i++)
        {
            for (j = i + 1; j < A->c; j++)
            {
                gr_swap(GR_MAT_ENTRY(A, i, j, sz), GR_MAT_ENTRY(A, j, i, sz), ctx);
            }
        }
    }
    else  /* Not aliased; general case */
    {
        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                status |= gr_set(GR_MAT_ENTRY(B, i, j, sz), GR_MAT_ENTRY(A, j, i, sz), ctx);
    }

    return status;
}
