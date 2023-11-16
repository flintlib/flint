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
_gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (res == poly)
    {
        slong i;

        for (i = 0; i < n / 2; i++)
            gr_swap(GR_ENTRY(res, i, sz), GR_ENTRY(res, n - 1 - i, sz), ctx);

        for (i = 0; i < n - len; i++)
            status |= gr_zero(GR_ENTRY(res, i, sz), ctx);
    }
    else
    {
        slong i;

        for (i = 0; i < n - len; i++)
            status |= gr_zero(GR_ENTRY(res, i, sz), ctx);

        for (i = 0; i < len; i++)
            status |= gr_set(GR_ENTRY(res, (n - len) + i, sz), GR_ENTRY(poly, (len - 1) - i, sz), ctx);
    }

    return status;
}

int
gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)
{
    slong len = FLINT_MIN(n, poly->length);
    int status;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    gr_poly_fit_length(res, n, ctx);
    status = _gr_poly_reverse(res->coeffs, poly->coeffs, len, n, ctx);
    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);

    return status;
}
