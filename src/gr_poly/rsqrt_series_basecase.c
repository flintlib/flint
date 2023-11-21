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

/* todo: use _gr_poly_pow_series_fmpq_recurrence if possible, when flen is short */
int
_gr_poly_rsqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (flen == 1)
    {
        status |= gr_rsqrt(res, f, ctx);
        if (status != GR_SUCCESS)
            return status;

        status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), len - 1, ctx);
    }
    else if (len == 2)
    {
        status |= gr_rsqrt(res, f, ctx);
        if (status != GR_SUCCESS)
            return status;

        status |= gr_mul(GR_ENTRY(res, 1, sz), res, GR_ENTRY(f, 1, sz), ctx);
        status |= gr_div(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), GR_ENTRY(f, 0, sz), ctx);
        status |= gr_mul_2exp_si(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), -1, ctx);
        status |= gr_neg(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), ctx);
    }
    else
    {
        gr_ptr t;
        GR_TMP_INIT_VEC(t, len, ctx);

        /* todo: when is it better and valid to invert first? */
        /* todo: somehow use gr_rsqrt() for the constant term? */
        /* todo: might want to call inv_series rather than inv_series_basecase here;
                 the cutoff could be much lower than that for sqrt */
        status |= _gr_poly_sqrt_series_basecase(t, f, flen, len, ctx);
        status |= _gr_poly_inv_series_basecase(res, t, len, len, ctx);

        GR_TMP_CLEAR_VEC(t, len, ctx);
    }

    return status;
}

int
gr_poly_rsqrt_series_basecase(gr_poly_t res, const gr_poly_t h, slong len, gr_ctx_t ctx)
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
        status = gr_poly_rsqrt_series_basecase(t, h, len, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_rsqrt_series_basecase(res->coeffs, h->coeffs, h->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
