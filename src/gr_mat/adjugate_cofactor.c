/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_adjugate_cofactor(gr_mat_t adj, gr_ptr det, const gr_mat_t A, gr_ctx_t ctx)
{
    gr_mat_t T;
    slong i, j, n, a, b;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx))
        return GR_DOMAIN;

    if (n == 0)
    {
        status |= gr_one(det, ctx);
        return status;
    }

    if (n == 1)
    {
        status |= gr_set(det, GR_MAT_ENTRY(A, 0, 0, sz), ctx);
        status |= gr_one(GR_MAT_ENTRY(adj, 0, 0, sz), ctx);
        return status;
    }

    if (n == 2)
    {
        gr_ptr t, u;
        GR_TMP_INIT2(t, u, ctx);
        status |= gr_mul(t, GR_MAT_ENTRY(A, 0, 0, sz), GR_MAT_ENTRY(A, 1, 1, sz), ctx);
        status |= gr_mul(u, GR_MAT_ENTRY(A, 0, 1, sz), GR_MAT_ENTRY(A, 1, 0, sz), ctx);
        status |= gr_set(GR_MAT_ENTRY(adj, 0, 0, sz), GR_MAT_ENTRY(A, 1, 1, sz), ctx);
        status |= gr_neg(GR_MAT_ENTRY(adj, 0, 1, sz), GR_MAT_ENTRY(A, 0, 1, sz), ctx);
        status |= gr_neg(GR_MAT_ENTRY(adj, 1, 0, sz), GR_MAT_ENTRY(A, 1, 0, sz), ctx);
        status |= gr_set(GR_MAT_ENTRY(adj, 1, 1, sz), GR_MAT_ENTRY(A, 0, 0, sz), ctx);
        status |= gr_sub(det, t, u, ctx);
        GR_TMP_CLEAR2(t, u, ctx);
        return status;
    }

    if (adj == A)
    {
        gr_mat_init(T, n, n, ctx);
        status |= gr_mat_adjugate_cofactor(T, det, A, ctx);
        gr_mat_swap(adj, T, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    /* todo: shallow copies */
    gr_mat_init(T, n - 1, n - 1, ctx);

    status |= gr_zero(det, ctx);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (a = 0; a < n; a++)
            {
                for (b = 0; b < n; b++)
                {
                    if (a != i && b != j)
                    {
                        status |= gr_set(GR_MAT_ENTRY(T, a - (a > i), b - (b > j), sz), GR_MAT_ENTRY(A, a, b, sz), ctx);
                    }
                }
            }

            status |= gr_mat_det(GR_MAT_ENTRY(adj, i, j, sz), T, ctx);
            if ((i + j) & 1)
                status |= gr_neg(GR_MAT_ENTRY(adj, i, j, sz), GR_MAT_ENTRY(adj, i, j, sz), ctx);

            if (i == 0)
            {
                status |= gr_addmul(det, GR_MAT_ENTRY(adj, i, j, sz), GR_MAT_ENTRY(A, i, j, sz), ctx);
            }
        }
    }

    status |= gr_mat_transpose(adj, adj, ctx);
    gr_mat_clear(T, ctx);

    return status;
}
