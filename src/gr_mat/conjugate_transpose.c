/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"

int
gr_mat_conjugate_transpose(gr_mat_t B, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(A, ctx);
    c = gr_mat_ncols(A, ctx);

    if (gr_mat_nrows(B, ctx) != c || gr_mat_ncols(B, ctx) != r)
        return GR_DOMAIN;

    if (A == B)
    {
        if (r == c)
        {
            for (i = 0; i < r; i++)
            {
                status |= gr_conj(GR_MAT_ENTRY(B, i, i, sz), GR_MAT_ENTRY(A, i, i, sz), ctx);

                for (j = i + 1; j < c; j++)
                {
                    status |= gr_conj(GR_MAT_ENTRY(B, i, j, sz), GR_MAT_ENTRY(A, i, j, sz), ctx);
                    status |= gr_conj(GR_MAT_ENTRY(B, j, i, sz), GR_MAT_ENTRY(A, j, i, sz), ctx);
                    gr_swap(GR_MAT_ENTRY(B, i, j, sz), GR_MAT_ENTRY(B, j, i, sz), ctx);
                }
            }
        }
        else
        {
            gr_mat_t T;
            gr_mat_init(T, c, r, ctx);
            status |= gr_mat_conjugate_transpose(T, A, ctx);
            gr_mat_swap(B, T, ctx);
            gr_mat_clear(T, ctx);
        }
    }
    else
    {
        for (i = 0; i < c; i++)
            for (j = 0; j < r; j++)
                status |= gr_conj(GR_MAT_ENTRY(B, i, j, sz), GR_MAT_ENTRY(A, j, i, sz), ctx);
    }

    return status;
}
