/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"

/* todo: allow overloading the following methods (or at least use
         vector functions)
   todo: quick bound versions */

int
gr_mat_norm_max(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr t;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    if (R == 0 || C == 0)
        return gr_zero(res, ctx);

    GR_TMP_INIT(t, ctx);

    for (i = 0; i < R; i++)
    {
        for (j = 0; j < C; j++)
        {
            if (i == 0 && j == 0)
                status |= gr_abs(res, GR_MAT_ENTRY(mat, i, j, sz), ctx);
            else
            {
                status |= gr_abs(t, GR_MAT_ENTRY(mat, i, j, sz), ctx);
                status |= gr_max(res, res, t, ctx);
            }
        }
    }

    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_mat_norm_1(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr s, t;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    if (R == 0 || C == 0)
        return gr_zero(res, ctx);

    GR_TMP_INIT2(s, t, ctx);

    for (j = 0; j < C; j++)
    {
        for (i = 0; i < R; i++)
        {
            if (i == 0)
                status |= gr_abs(s, GR_MAT_ENTRY(mat, i, j, sz), ctx);
            else
            {
                status |= gr_abs(t, GR_MAT_ENTRY(mat, i, j, sz), ctx);
                status |= gr_add(s, s, t, ctx);
            }
        }

        if (j == 0)
            gr_swap(res, s, ctx);
        else
            status |= gr_max(res, res, s, ctx);
    }

    GR_TMP_CLEAR2(s, t, ctx);

    return status;
}

int
gr_mat_norm_inf(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr s, t;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    if (R == 0 || C == 0)
        return gr_zero(res, ctx);

    GR_TMP_INIT2(s, t, ctx);

    for (i = 0; i < R; i++)
    {
        for (j = 0; j < C; j++)
        {
            if (j == 0)
                status |= gr_abs(s, GR_MAT_ENTRY(mat, i, j, sz), ctx);
            else
            {
                status |= gr_abs(t, GR_MAT_ENTRY(mat, i, j, sz), ctx);
                status |= gr_add(s, s, t, ctx);
            }
        }

        if (i == 0)
            gr_swap(res, s, ctx);
        else
            status |= gr_max(res, res, s, ctx);
    }

    GR_TMP_CLEAR2(s, t, ctx);

    return status;
}

int
gr_mat_norm_frobenius(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr t;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    if (R == 0 || C == 0)
        return gr_zero(res, ctx);

    GR_TMP_INIT(t, ctx);

    status |= gr_zero(res, ctx);

    for (i = 0; i < R; i++)
    {
        for (j = 0; j < C; j++)
        {
            status |= gr_abs(t, GR_MAT_ENTRY(mat, i, j, sz), ctx);
            status |= gr_sqr(t, t, ctx);
            status |= gr_add(res, res, t, ctx);
        }
    }

    status |= gr_sqrt(res, res, ctx);

    GR_TMP_CLEAR(t, ctx);

    return status;
}
