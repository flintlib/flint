/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
_gr_poly_evaluate_other_horner(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len == 0)
    {
        return gr_zero(res, x_ctx);
    }
    else if (len == 1 || gr_is_zero(x, x_ctx) == T_TRUE)
    {
        return gr_set_other(res, f, ctx, x_ctx);
    }
    else if (len == 2)
    {
        slong sz = ctx->sizeof_elem;

        status |= gr_mul_other(res, x, GR_ENTRY(f, 1, sz), ctx, x_ctx);
        status |= gr_add_other(res, res, f, ctx, x_ctx);

        return status;
    }
    else
    {
        slong i = len - 1;
        slong sz = ctx->sizeof_elem;
        gr_ptr t, u;

        GR_TMP_INIT2(t, u, x_ctx);

        status |= gr_set_other(u, GR_ENTRY(f, i, sz), ctx, x_ctx);

        for (i = len - 2; i >= 0; i--)
        {
            status |= gr_mul(t, u, x, x_ctx);
            status |= gr_add_other(u, t, GR_ENTRY(f, i, sz), ctx, x_ctx);
        }

        gr_swap(res, u, x_ctx);

        GR_TMP_CLEAR2(t, u, x_ctx);

        return status;
    }
}

int
gr_poly_evaluate_other_horner(gr_ptr res, const gr_poly_t f, gr_srcptr a, gr_ctx_t a_ctx, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_other_horner(res, f->coeffs, f->length, a, a_ctx, ctx);
}
