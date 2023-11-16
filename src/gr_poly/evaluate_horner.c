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
_gr_poly_evaluate_horner(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len == 0)
    {
        return gr_zero(res, ctx);
    }
    else if (len == 1 || gr_is_zero(x, ctx) == T_TRUE)
    {
        return gr_set(res, f, ctx);
    }
    else if (len == 2)
    {
        slong sz = ctx->sizeof_elem;

        status |= gr_mul(res, x, GR_ENTRY(f, 1, sz), ctx);
        status |= gr_add(res, res, f, ctx);

        return status;
    }
    else
    {
        slong i = len - 1;
        slong sz = ctx->sizeof_elem;
        gr_ptr t, u;

        GR_TMP_INIT2(t, u, ctx);

        status |= gr_set(u, GR_ENTRY(f, i, sz), ctx);

        for (i = len - 2; i >= 0; i--)
        {
            status |= gr_mul(t, u, x, ctx);
            status |= gr_add(u, GR_ENTRY(f, i, sz), t, ctx);
        }

        gr_swap(res, u, ctx);

        GR_TMP_CLEAR2(t, u, ctx);

        return status;
    }
}

int
gr_poly_evaluate_horner(gr_ptr res, const gr_poly_t f, gr_srcptr a, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_horner(res, f->coeffs, f->length, a, ctx);
}
