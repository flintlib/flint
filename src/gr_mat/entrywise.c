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

int
gr_mat_entrywise_unary_op(gr_mat_t res, gr_method_unary_op f, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    if (R != gr_mat_nrows(res, ctx) || C != gr_mat_ncols(res, ctx))
        return GR_DOMAIN;

    for (i = 0; i < R; i++)
        for (j = 0; j < C; j++)
            status |= f(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);

    return status;
}

int
gr_mat_entrywise_binary_op(gr_mat_t res, gr_method_binary_op f, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    R = gr_mat_nrows(mat1, ctx);
    C = gr_mat_ncols(mat1, ctx);

    if (R != gr_mat_nrows(res, ctx) || C != gr_mat_ncols(res, ctx) || R != gr_mat_nrows(mat2, ctx) || C != gr_mat_ncols(mat2, ctx))
        return GR_DOMAIN;

    for (i = 0; i < R; i++)
        for (j = 0; j < C; j++)
            status |= f(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat1, i, j, sz), GR_MAT_ENTRY(mat2, i, j, sz), ctx);

    return status;
}

int
gr_mat_entrywise_binary_op_scalar(gr_mat_t res, gr_method_binary_op f, const gr_mat_t mat, gr_srcptr c, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    if (R != gr_mat_nrows(res, ctx) || C != gr_mat_ncols(res, ctx))
        return GR_DOMAIN;

    for (i = 0; i < R; i++)
        for (j = 0; j < C; j++)
            status |= f(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), c, ctx);

    return status;
}

truth_t
gr_mat_entrywise_unary_predicate_all(gr_method_unary_predicate f, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    truth_t val, ans = T_TRUE;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    for (i = 0; i < R; i++)
    {
        for (j = 0; j < C; j++)
        {
            val = f(GR_MAT_ENTRY(mat, i, j, sz), ctx);
            if (val == T_FALSE)
                return T_FALSE;
            ans = truth_and(ans, val);
        }
    }

    return ans;
}

truth_t
gr_mat_entrywise_unary_predicate_any(gr_method_unary_predicate f, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    truth_t val, ans = T_FALSE;

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    for (i = 0; i < R; i++)
    {
        for (j = 0; j < C; j++)
        {
            val = f(GR_MAT_ENTRY(mat, i, j, sz), ctx);
            if (val == T_TRUE)
                return T_TRUE;
            ans = truth_or(ans, val);
        }
    }

    return ans;
}

truth_t
gr_mat_entrywise_binary_predicate_all(gr_method_binary_predicate f, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    slong R, C, i, j;
    slong sz = ctx->sizeof_elem;
    truth_t val, ans = T_TRUE;

    R = gr_mat_nrows(mat1, ctx);
    C = gr_mat_ncols(mat1, ctx);

    if (R != gr_mat_nrows(mat2, ctx) || C != gr_mat_ncols(mat2, ctx))
        return T_FALSE;

    for (i = 0; i < R; i++)
    {
        for (j = 0; j < C; j++)
        {
            val = f(GR_MAT_ENTRY(mat1, i, j, sz), GR_MAT_ENTRY(mat2, i, j, sz), ctx);
            if (val == T_FALSE)
                return T_FALSE;
            ans = truth_and(ans, val);
        }
    }

    return ans;
}
