/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

/* todo: overloadable function */
int
_gr_poly_derivative(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op_ui mul_ui = GR_BINARY_OP_UI(ctx, MUL_UI);
    slong i;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    for (i = 1; i < len; i++)
        status |= mul_ui(GR_ENTRY(res, i - 1, sz), GR_ENTRY(poly, i, sz), i, ctx);

    return status;
}

int
gr_poly_derivative(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len <= 1)
    {
        status = gr_poly_zero(res, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len - 1, ctx);
        status = _gr_poly_derivative(res->coeffs, poly->coeffs, len, ctx);
        _gr_poly_set_length(res, len - 1, ctx);
        /* todo: only call in nonzero characteristic */
        _gr_poly_normalise(res, ctx);
    }

    return status;
}
