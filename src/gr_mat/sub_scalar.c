/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"

int
gr_mat_sub_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (res == mat)
    {
        for (i = 0; i < FLINT_MIN(r, c); i++)
            status |= gr_sub(GR_MAT_ENTRY(res, i, i, sz), GR_MAT_ENTRY(res, i, i, sz), x, ctx);
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            for (j = 0; j < c; j++)
            {
                /* todo: vectorize */
                if (i == j)
                    status |= gr_sub(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), x, ctx);
                else
                    status |= gr_set(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);
            }
        }
    }

    return status;
}

int
gr_mat_scalar_sub(gr_mat_t res, gr_srcptr x, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            /* todo: vectorize */
            if (i == j)
                status |= gr_sub(GR_MAT_ENTRY(res, i, j, sz), x, GR_MAT_ENTRY(mat, i, j, sz), ctx);
            else
                status |= gr_neg(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);
        }
    }

    return status;
}

int
gr_mat_sub_scalar_other(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (res == mat)
    {
        for (i = 0; i < FLINT_MIN(r, c); i++)
            status |= gr_sub_other(GR_MAT_ENTRY(res, i, i, sz), GR_MAT_ENTRY(res, i, i, sz), x, x_ctx, ctx);
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            for (j = 0; j < c; j++)
            {
                /* todo: vectorize */
                if (i == j)
                    status |= gr_sub_other(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), x, x_ctx, ctx);
                else
                    status |= gr_set(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);
            }
        }
    }

    return status;
}

int
gr_mat_scalar_other_sub(gr_mat_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            /* todo: vectorize */
            if (i == j)
                status |= gr_other_sub(GR_MAT_ENTRY(res, i, j, sz), x, x_ctx, GR_MAT_ENTRY(mat, i, j, sz), ctx);
            else
                status |= gr_neg(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);
        }
    }

    return status;
}

int
gr_mat_sub_ui(gr_mat_t res, const gr_mat_t mat, ulong x, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (res == mat)
    {
        for (i = 0; i < FLINT_MIN(r, c); i++)
            status |= gr_sub_ui(GR_MAT_ENTRY(res, i, i, sz), GR_MAT_ENTRY(res, i, i, sz), x, ctx);
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            for (j = 0; j < c; j++)
            {
                /* todo: vectorize */
                if (i == j)
                    status |= gr_sub_ui(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), x, ctx);
                else
                    status |= gr_set(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);
            }
        }
    }

    return status;
}

int
gr_mat_sub_si(gr_mat_t res, const gr_mat_t mat, slong x, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (res == mat)
    {
        for (i = 0; i < FLINT_MIN(r, c); i++)
            status |= gr_sub_si(GR_MAT_ENTRY(res, i, i, sz), GR_MAT_ENTRY(res, i, i, sz), x, ctx);
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            for (j = 0; j < c; j++)
            {
                /* todo: vectorize */
                if (i == j)
                    status |= gr_sub_si(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), x, ctx);
                else
                    status |= gr_set(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);
            }
        }
    }

    return status;
}

int
gr_mat_sub_fmpz(gr_mat_t res, const gr_mat_t mat, const fmpz_t x, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (res == mat)
    {
        for (i = 0; i < FLINT_MIN(r, c); i++)
            status |= gr_sub_fmpz(GR_MAT_ENTRY(res, i, i, sz), GR_MAT_ENTRY(res, i, i, sz), x, ctx);
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            for (j = 0; j < c; j++)
            {
                /* todo: vectorize */
                if (i == j)
                    status |= gr_sub_fmpz(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), x, ctx);
                else
                    status |= gr_set(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);
            }
        }
    }

    return status;
}

int
gr_mat_sub_fmpq(gr_mat_t res, const gr_mat_t mat, const fmpq_t x, gr_ctx_t ctx)
{
    slong i, j, r, c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (res == mat)
    {
        for (i = 0; i < FLINT_MIN(r, c); i++)
            status |= gr_sub_fmpq(GR_MAT_ENTRY(res, i, i, sz), GR_MAT_ENTRY(res, i, i, sz), x, ctx);
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            for (j = 0; j < c; j++)
            {
                /* todo: vectorize */
                if (i == j)
                    status |= gr_sub_fmpq(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), x, ctx);
                else
                    status |= gr_set(GR_MAT_ENTRY(res, i, j, sz), GR_MAT_ENTRY(mat, i, j, sz), ctx);
            }
        }
    }

    return status;
}
