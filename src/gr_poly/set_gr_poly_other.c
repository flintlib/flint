/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
gr_poly_set_gr_poly_other(gr_poly_t res, const gr_poly_t x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong x_sz = x_ctx->sizeof_elem;
    slong sz = ctx->sizeof_elem;
    slong i, len = gr_poly_length(x, x_ctx);

    gr_poly_fit_length(res, len, ctx);
    _gr_poly_set_length(res, len, ctx);

    for (i = 0; i < len; i++)
    {
        status |= gr_set_other(GR_ENTRY(res->coeffs, i, sz), GR_ENTRY(x->coeffs, i, x_sz), x_ctx, ctx);
    }

    if (status == GR_SUCCESS)
        _gr_poly_normalise(res, ctx);
    else
        _gr_poly_set_length(res, 0, ctx);

    return status;
}
