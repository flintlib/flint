/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
_gr_poly_integral(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    slong k = len - 1;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    for (k = len - 1; k > 0 && status == GR_SUCCESS; k--)
        status |= gr_div_ui(GR_ENTRY(res, k, sz), GR_ENTRY(poly, k - 1, sz), k, ctx);

    status |= gr_zero(res, ctx);
    return status;
}

int
gr_poly_integral(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)
{
    int status;
    gr_poly_fit_length(res, poly->length + 1, ctx);
    status = _gr_poly_integral(res->coeffs, poly->coeffs, poly->length + 1, ctx);
    _gr_poly_set_length(res, poly->length + 1, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
