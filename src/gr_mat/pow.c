/*
    Copyright (C) 2025 Lars Göttgens
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fmpz.h"
#include "fmpq.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_mat.h"

int
gr_mat_pow_ui(gr_mat_t res, const gr_mat_t mat, ulong exp, gr_ctx_t ctx)
{
    int status;
    slong sz = ctx->sizeof_elem;
    slong d;

    d = gr_mat_nrows(res, ctx);

    if (d != gr_mat_ncols(res, ctx) || d != gr_mat_nrows(mat, ctx)
        || d != gr_mat_ncols(mat, ctx))
    {
        return GR_DOMAIN;
    }

    status = GR_SUCCESS;

    if (exp <= 2 || d <= 1)
    {
        if (exp == 0 || d == 0)
        {
            status |= gr_mat_one(res, ctx);
        }
        else if (d == 1)
        {
            status |= gr_pow_ui(GR_MAT_ENTRY(res, 0, 0, sz),
                                GR_MAT_ENTRY(mat, 0, 0, sz), exp, ctx);
        }
        else if (exp == 1)
        {
            status |= gr_mat_set(res, mat, ctx);
        }
        else if (exp == 2)
        {
            status |= gr_mat_sqr(res, mat, ctx);
        }
    }
    else
    {
        gr_ctx_t mctx;

        gr_ctx_init_matrix_ring(mctx, ctx, d);
        status |= gr_generic_pow_ui(res, mat, exp, mctx);
        gr_ctx_clear(mctx);
    }

    return status;
}

int
gr_mat_pow_fmpz(gr_mat_t res, const gr_mat_t mat, const fmpz_t exp, gr_ctx_t ctx)
{
    int status;
    slong sz = ctx->sizeof_elem;
    slong d;

    d = gr_mat_nrows(res, ctx);

    if (d != gr_mat_ncols(res, ctx) || d != gr_mat_nrows(mat, ctx)
        || d != gr_mat_ncols(mat, ctx))
    {
        return GR_DOMAIN;
    }

    status = GR_SUCCESS;


    if (fmpz_is_zero(exp) || d == 0)
    {
        status |= gr_mat_one(res, ctx);
    }
    else if (d == 1)
    {
        status |= gr_pow_fmpz(GR_MAT_ENTRY(res, 0, 0, sz),
                              GR_MAT_ENTRY(mat, 0, 0, sz), exp, ctx);
    }
    else if (fmpz_is_one(exp))
    {
        status |= gr_mat_set(res, mat, ctx);
    }
    else
    {
        gr_ctx_t mctx;

        gr_ctx_init_matrix_ring(mctx, ctx, d);
        status |= gr_generic_pow_fmpz(res, mat, exp, mctx);
        gr_ctx_clear(mctx);
    }

    return status;
}

int
gr_mat_pow_si(gr_mat_t res, const gr_mat_t mat, slong exp, gr_ctx_t ctx)
{
    if (exp >= 0)
    {
        return gr_mat_pow_ui(res, mat, exp, ctx);
    }
    else
    {
        fmpz_t n;
        int status;
        fmpz_init_set_si(n, exp);
        status = gr_mat_pow_fmpz(res, mat, n, ctx);
        fmpz_clear(n);
        return status;
    }
}

int
gr_pow_jet(gr_ptr res, gr_srcptr x, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    slong i;
    int status;
    slong sz = ctx->sizeof_elem;

    if (len <= 0)
        return GR_SUCCESS;

    status = gr_pow(res, x, c, ctx);

    if (status == GR_SUCCESS && len > 1)
    {
        gr_ptr t, u;
        GR_TMP_INIT2(t, u, ctx);

        status |= gr_inv(t, x, ctx);

        for (i = 1; i < len; i++)
        {
            status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), t, ctx);
            status |= gr_sub_ui(u, c, i - 1, ctx);
            status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), u, ctx);
            status |= gr_div_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), i, ctx);
        }

        GR_TMP_CLEAR(t, ctx);
        GR_TMP_CLEAR(u, ctx);
    }

    return status;
}

int
gr_mat_pow_scalar_jordan(gr_mat_t res, const gr_mat_t A, gr_srcptr c, gr_ctx_t ctx)
{
    return gr_mat_func_param_jordan(res, A, (gr_method_vec_scalar_op) gr_pow_jet, c, ctx);
}

int
gr_mat_pow_scalar(gr_mat_t res, const gr_mat_t A, gr_srcptr c, gr_ctx_t ctx)
{
    slong n;

    /* we don't look for fmpz because a floating-point exponent
       could be something huge */
    if (gr_get_si(&n, c, ctx) == GR_SUCCESS)
    {
        return gr_mat_pow_si(res, A, n, ctx);
    }
    else
    {
        int status = gr_mat_pow_scalar_jordan(res, A, c, ctx);

        /* We cannot conclude nonexistence of a power just
           because the Jordan algorithm failed. */
        return (status == GR_SUCCESS) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
gr_pow_fmpq_jet(gr_ptr res, gr_srcptr x, slong len, const fmpq_t c, gr_ctx_t ctx)
{
    slong i;
    int status;
    slong sz = ctx->sizeof_elem;

    if (len <= 0)
        return GR_SUCCESS;

    status = gr_pow_fmpq(res, x, c, ctx);

    if (status == GR_SUCCESS && len > 1)
    {
        gr_ptr t;
        fmpq_t u;

        GR_TMP_INIT(t, ctx);
        fmpq_init(u);

        status |= gr_inv(t, x, ctx);

        for (i = 1; i < len; i++)
        {
            status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), t, ctx);
            fmpq_sub_ui(u, c, i - 1);
            status |= gr_mul_fmpq(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), u, ctx);
            status |= gr_div_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), i, ctx);
        }

        GR_TMP_CLEAR(t, ctx);
        fmpq_clear(u);
    }

    return status;
}

int
gr_mat_pow_fmpq_jordan(gr_mat_t res, const gr_mat_t A, const fmpq_t c, gr_ctx_t ctx)
{
    return gr_mat_func_param_jordan(res, A, (gr_method_vec_scalar_op) gr_pow_fmpq_jet, c, ctx);
}

int
gr_mat_pow_fmpq(gr_mat_t res, const gr_mat_t A, const fmpq_t c, gr_ctx_t ctx)
{
    if (fmpz_is_one(fmpq_denref(c)))
        return gr_mat_pow_fmpz(res, A, fmpq_numref(c), ctx);
    else
    {
        int status = gr_mat_pow_fmpq_jordan(res, A, c, ctx);

        /* We cannot conclude nonexistence of an Nth root just
           because the Jordan algorithm failed. */
        return (status == GR_SUCCESS) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
gr_mat_sqrt(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    fmpq_t c;
    *fmpq_numref(c) = 1;
    *fmpq_denref(c) = 2;
    return gr_mat_pow_fmpq(res, A, c, ctx);
}

int
gr_mat_rsqrt(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    fmpq_t c;
    *fmpq_numref(c) = -1;
    *fmpq_denref(c) = 2;
    return gr_mat_pow_fmpq(res, A, c, ctx);
}
