/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: Karp-Markstein */
int
_gr_poly_sqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ptr t;

    GR_TMP_INIT_VEC(t, len, ctx);

    status = _gr_poly_rsqrt_series_newton(t, f, flen, len, cutoff, ctx);
    if (status == GR_SUCCESS)
        status |= _gr_poly_mullow(res, t, len, f, flen, len, ctx);

    GR_TMP_CLEAR_VEC(t, len, ctx);

    return status;
}

int
gr_poly_sqrt_series_newton(gr_poly_t res, const gr_poly_t h, slong len, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen;

    hlen = h->length;

    if (len == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (hlen == 1)
        len = 1;

    if (res == h)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_sqrt_series_newton(t, h, len, cutoff, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_sqrt_series_newton(res->coeffs, h->coeffs, h->length, len, cutoff, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
