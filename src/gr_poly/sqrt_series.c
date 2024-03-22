/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

/* todo: characteristic 2, ... */

static truth_t
_char_two(gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t res;
    GR_TMP_INIT(t, ctx);
    GR_MUST_SUCCEED(gr_set_ui(t, 2, ctx));
    res = gr_is_zero(t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return res;
}

int
_gr_poly_sqrt_series_generic(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    int status;

    status = _gr_poly_sqrt_series_newton(res, f, flen, len, 2, ctx);

    if (status == GR_DOMAIN && (gr_ctx_is_field(ctx) != T_TRUE || _char_two(ctx) != T_FALSE))
        return GR_UNABLE;

    return status;
}

int
gr_poly_sqrt_series(gr_poly_t res, const gr_poly_t h, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen;

    hlen = h->length;

    if (hlen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (hlen == 1)
        len = 1;

    if (res == h)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_sqrt_series(t, h, len, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_sqrt_series(res->coeffs, h->coeffs, h->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
