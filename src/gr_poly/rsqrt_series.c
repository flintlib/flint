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
_gr_poly_rsqrt_series_generic(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    if (flen <= 8 || ctx->methods[GR_METHOD_POLY_MULLOW] == (gr_funcptr) _gr_poly_mullow_generic)
        return _gr_poly_rsqrt_series_basecase(res, f, flen, len, ctx);
    else
        return _gr_poly_rsqrt_series_newton(res, f, flen, len, FLINT_MIN(10, len / 2), ctx);
}

int
gr_poly_rsqrt_series(gr_poly_t res, const gr_poly_t h, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    hlen = h->length;

    if (hlen == 0)
        return GR_DOMAIN;

    if (hlen == 1)
        len = 1;

    if (res == h)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_rsqrt_series(t, h, len, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_rsqrt_series(res->coeffs, h->coeffs, h->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
