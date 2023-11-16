/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

/* todo: move me */
int
gr_dot_other(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx2, gr_ctx_t ctx)
{
    int status;
    slong i, sz, sz2;
    gr_ptr t;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    sz2 = ctx2->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    if (initial == NULL)
    {
        status |= gr_mul_other(res, vec1, vec2, ctx2, ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= gr_mul_other(t, vec1, vec2, ctx2, ctx);
        status |= gr_add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= gr_mul_other(t, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz2), ctx2, ctx);
        status |= gr_add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int
_gr_poly_evaluate_other_rectangular(gr_ptr y, gr_srcptr poly,
    slong len, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong x_sz = x_ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len <= 2)
    {
        if (len == 0)
            return gr_zero(y, x_ctx);

        if (len == 1)
            return gr_set_other(y, poly, ctx, x_ctx);

        status |= gr_mul_other(y, x, GR_ENTRY(poly, 1, sz), ctx, x_ctx);
        status |= gr_add_other(y, y, GR_ENTRY(poly, 0, sz), ctx, x_ctx);
        return status;
    }
    else
    {
        slong i, m, r;
        gr_ptr xs;
        gr_ptr s, t, c;

        m = n_sqrt(len) + 1;
        r = (len + m - 1) / m;

        GR_TMP_INIT_VEC(xs, m + 1, x_ctx);
        GR_TMP_INIT3(s, t, c, x_ctx);

        status |= _gr_vec_set_powers(xs, x, m + 1, x_ctx);
        status |= gr_dot_other(y, NULL, 0, GR_ENTRY(xs, 1, x_sz), GR_ENTRY(poly, (r - 1) * m + 1, sz), len - (r - 1) * m - 1, ctx, x_ctx);
        status |= gr_add_other(y, y, GR_ENTRY(poly, (r - 1) * m, sz), ctx, x_ctx);

        for (i = r - 2; i >= 0; i--)
        {
            status |= gr_dot_other(s, NULL, 0, GR_ENTRY(xs, 1, x_sz), GR_ENTRY(poly, i * m + 1, sz), m - 1, ctx, x_ctx);
            status |= gr_add_other(s, s, GR_ENTRY(poly, i * m, sz), ctx, x_ctx);
            status |= gr_mul(y, y, GR_ENTRY(xs, m, x_sz), x_ctx);
            status |= gr_add(y, y, s, x_ctx);
        }

        GR_TMP_CLEAR_VEC(xs, m + 1, x_ctx);
        GR_TMP_CLEAR3(s, t, c, x_ctx);

        return status;
    }
}

int
gr_poly_evaluate_other_rectangular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_other_rectangular(res, f->coeffs, f->length, x, x_ctx, ctx);
}
