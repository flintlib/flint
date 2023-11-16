/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: generalize */
#include "arb.h"
#include "acb.h"

static int
want_division(gr_srcptr f, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_RR_ARB)
        return arb_bits(f) <= 0.25 * _gr_ctx_get_real_prec(ctx);

    if (ctx->which_ring == GR_CTX_CC_ACB)
        return acb_bits(f) <= 0.25 * _gr_ctx_get_real_prec(ctx);

    return 0;
}

int
_gr_poly_taylor_shift_convolution(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, n = len - 1;
    gr_ptr f, d, t, u;

    if (res != poly)
        status |= _gr_vec_set(res, poly, len, ctx);

    if (gr_is_zero(c, ctx) == T_TRUE  || len <= 1)
        return status;

    GR_TMP_INIT_VEC(t, 2 * len + 2, ctx);
    u = GR_ENTRY(t, len, sz);
    f = GR_ENTRY(u, len, sz);
    d = GR_ENTRY(f, 1, sz);

    status |= gr_one(f, ctx);

    /* todo: vector function for Borel transform */
    for (i = 2; i <= n; i++)
    {
        status |= gr_mul_ui(f, f, i, ctx);
        status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), f, ctx);
    }

    status |= _gr_poly_reverse(res, res, len, len, ctx);

    status |= gr_one(GR_ENTRY(t, n, sz), ctx);
    for (i = n; i > 0; i--)
        status |= gr_mul_ui(GR_ENTRY(t, i - 1, sz), GR_ENTRY(t, i, sz), i, ctx);

    if (gr_is_neg_one(c, ctx) == T_TRUE)
    {
        for (i = 1; i <= n; i += 2)
            status |= gr_neg(GR_ENTRY(t, i, sz), GR_ENTRY(t, i, sz), ctx);
    }
    else if (gr_is_one(c, ctx) != T_TRUE)
    {
        status |= gr_set(d, c, ctx);

        for (i = 1; i <= n; i++)
        {
            status |= gr_mul(GR_ENTRY(t, i, sz), GR_ENTRY(t, i, sz), d, ctx);
            status |= gr_mul(d, d, c, ctx);
        }
    }

    status |= _gr_poly_mullow(u, res, len, t, len, len, ctx);
    status |= gr_mul(f, f, f, ctx);

    if (want_division(f, ctx))
    {
        for (i = 0; i <= n; i++)
            status |= gr_div(GR_ENTRY(u, i, sz), GR_ENTRY(u, i, sz), f, ctx);

        status |= gr_one(f, ctx);
    }
    else
    {
        status |= gr_inv(f, f, ctx);
    }

    for (i = n; i >= 0; i--)
    {
        status |= gr_mul(GR_ENTRY(res, i, sz), GR_ENTRY(u, n - i, sz), f, ctx);
        status |= gr_mul_ui(f, f, (i == 0) ? 1 : i, ctx);
    }

    GR_TMP_CLEAR_VEC(t, 2 * len + 2, ctx);

    return status;
}

int
gr_poly_taylor_shift_convolution(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (res != f)
        status |= gr_poly_set(res, f, ctx);

    status |= _gr_poly_taylor_shift_convolution(res->coeffs, res->coeffs, res->length, c, ctx);
    return status;
}
