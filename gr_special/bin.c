/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_special.h"

/* todo: other algorithms; specializations */

int
gr_rising_ui_forward(gr_ptr res, gr_srcptr x, ulong n, gr_ctx_t ctx)
{
    gr_ptr t;
    ulong k;
    int status = GR_SUCCESS;

    if (n <= 1)
    {
        if (n == 0)
            return gr_one(res, ctx);
        else
            return gr_set(res, x, ctx);
    }

    GR_TMP_INIT(t, ctx);

    status |= gr_add_ui(t, x, 1, ctx);
    status |= gr_mul(res, x, t, ctx);

    for (k = 2; k < n; k++)
    {
        status |= gr_add_ui(t, t, 1, ctx);
        status |= gr_mul(res, res, t, ctx);
    }

    GR_TMP_CLEAR(t, ctx);

    return status;
}

static int
bsplit(gr_ptr y, gr_srcptr x, ulong a, ulong b, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (b - a <= 16)
    {
        if (a == 0)
        {
            status |= gr_rising_ui_forward(y, x, b, ctx);
        }
        else
        {
            status |= gr_add_ui(y, x, a, ctx);
            status |= gr_rising_ui_forward(y, y, b - a, ctx);
        }
    }
    else
    {
        gr_ptr t, u;
        ulong m = a + (b - a) / 2;

        GR_TMP_INIT2(t, u, ctx);

        status |= bsplit(t, x, a, m, ctx);
        status |= bsplit(u, x, m, b, ctx);
        status |= gr_mul(y, t, u, ctx);

        GR_TMP_CLEAR2(t, u, ctx);
    }

    return status;
}

int gr_generic_rising_ui(gr_ptr res, gr_srcptr x, ulong n, gr_ctx_t ctx)
{
    return bsplit(res, x, 0, n, ctx);
}

int gr_generic_falling_ui(gr_ptr res, gr_srcptr x, ulong n, gr_ctx_t ctx)
{
    if (n == 0)
        return gr_one(res, ctx);
    else
        return gr_sub_ui(res, x, n - 1, ctx) | gr_rising_ui(res, res, n, ctx);
}

int gr_generic_rising(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    ulong n;

    if (gr_get_ui(&n, y, ctx) == GR_SUCCESS)
        return gr_rising_ui(res, x, n, ctx);

    return GR_UNABLE;
}

int gr_generic_falling(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_sub_ui(t, y, 1, ctx);
    status |= gr_sub(t, x, t, ctx);
    status |= gr_rising(res, t, y, ctx);

    GR_TMP_CLEAR(t, ctx);

    return status;
}

/* todo: addition basecase */

int
gr_generic_bin_ui_vec(gr_ptr res, ulong n, slong len, gr_ctx_t ctx)
{
    gr_ptr f;
    slong sz = ctx->sizeof_elem;
    slong i, m;
    int status = GR_SUCCESS;
    truth_t finite_char;

    if (len <= 0)
        return GR_SUCCESS;

    if (len == 1)
        return gr_one(res, ctx);

    m = FLINT_MIN(len - 1, n / 2) + 1;

    finite_char = gr_ctx_is_finite_characteristic(ctx);

    if (finite_char == T_TRUE)
    {
        status = _gr_vec_reciprocals(GR_ENTRY(res, 1, sz), m - 1, ctx);

        if (status == GR_SUCCESS)
        {
            gr_method_binary_op_ui mul_ui = GR_BINARY_OP_UI(ctx, MUL_UI);
            gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);

            status |= gr_one(res, ctx);

            for (i = 1; i < m; i++)
            {
                status |= mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), ctx);
                status |= mul_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), n - i + 1, ctx);
            }
        }
        else
        {
            GR_TMP_INIT_VEC(f, 2, ctx);

            status = gr_one(f, ctx);
            status |= gr_one(GR_ENTRY(f, 1, sz), ctx);
            status |= _gr_poly_pow_series_ui_binexp(res, f, 2, n, m, ctx);

            GR_TMP_CLEAR_VEC(f, 2, ctx);
        }
    }
    else
    {
        gr_method_binary_op_ui mul_ui = GR_BINARY_OP_UI(ctx, MUL_UI);
        gr_method_binary_op_ui div_ui = GR_BINARY_OP_UI(ctx, DIV_UI);
        gr_method_binary_op_ui divexact_ui = GR_BINARY_OP_UI(ctx, DIVEXACT_UI);

        status |= gr_one(res, ctx);

        if (finite_char == T_FALSE)
        {
            for (i = 1; i < m; i++)
            {
                status |= mul_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), n - i + 1, ctx);
                status |= divexact_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), i, ctx);
            }
        }
        else
        {
            for (i = 1; i < m; i++)
            {
                status |= mul_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), n - i + 1, ctx);
                status |= div_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), i, ctx);
            }
        }
    }

    /* todo: vec reverse */
    if (m < len)
    {
        for (i = m; i <= FLINT_MIN(n, len - 1); i++)
            status |= gr_set(GR_ENTRY(res, i, sz), GR_ENTRY(res, n - i, sz), ctx);
    }

    if (len - 1 > n)
        status |= _gr_vec_zero(GR_ENTRY(res, n + 1, sz), len - 1 - n, ctx);

    return status;
}
