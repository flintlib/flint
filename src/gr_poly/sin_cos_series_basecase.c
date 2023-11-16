/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_sin_cos_series_basecase(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen,
        slong n, int times_pi, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong k, alen = FLINT_MIN(n, hlen);
    gr_ptr a, t, u;

    hlen = FLINT_MIN(hlen, n);

    if (times_pi)
        status |= gr_sin_cos_pi(s, c, h, ctx);
    else
        status |= gr_sin_cos(s, c, h, ctx);

    if (hlen == 1)
    {
        status |= _gr_vec_zero(GR_ENTRY(s, 1, sz), n - 1, ctx);
        status |= _gr_vec_zero(GR_ENTRY(c, 1, sz), n - 1, ctx);
        return status;
    }

    GR_TMP_INIT_VEC(a, alen + 2, ctx);
    t = GR_ENTRY(a, alen, sz);
    u = GR_ENTRY(t, 1, sz);

    for (k = 1; k < alen; k++)
        status |= gr_mul_ui(GR_ENTRY(a, k, sz), GR_ENTRY(h, k, sz), k, ctx);

    if (times_pi)
    {
        status |= gr_pi(t, ctx);
        status |= _gr_vec_mul_scalar(GR_ENTRY(a, 1, sz), GR_ENTRY(a, 1, sz), alen - 1, t, ctx);
    }

    for (k = 1; k < n; k++)
    {
        slong l = FLINT_MIN(k, hlen - 1);

        status |= _gr_vec_dot_rev(t, NULL, 1, GR_ENTRY(a, 1, sz), GR_ENTRY(s, k - l, sz), l, ctx);
        status |= _gr_vec_dot_rev(u, NULL, 0, GR_ENTRY(a, 1, sz), GR_ENTRY(c, k - l, sz), l, ctx);

        status |= gr_div_ui(GR_ENTRY(c, k, sz), t, k, ctx);
        status |= gr_div_ui(GR_ENTRY(s, k, sz), u, k, ctx);
    }

    GR_TMP_CLEAR_VEC(a, alen + 2, ctx);

    return status;
}

int
gr_poly_sin_cos_series_basecase(gr_poly_t s, gr_poly_t c,
        const gr_poly_t h, slong n, int times_pi, gr_ctx_t ctx)
{
    slong hlen = h->length;
    int status = GR_SUCCESS;

    if (n == 0)
    {
        status |= gr_poly_zero(s, ctx);
        status |= gr_poly_zero(c, ctx);
        return status;
    }

    if (hlen == 0)
    {
        status |= gr_poly_zero(s, ctx);
        status |= gr_poly_one(c, ctx);
        return status;
    }

    if (hlen == 1)
        n = 1;

    gr_poly_fit_length(s, n, ctx);
    gr_poly_fit_length(c, n, ctx);
    status |= _gr_poly_sin_cos_series_basecase(s->coeffs, c->coeffs, h->coeffs, hlen, n, times_pi, ctx);
    _gr_poly_set_length(s, n, ctx);
    _gr_poly_normalise(s, ctx);
    _gr_poly_set_length(c, n, ctx);
    _gr_poly_normalise(c, ctx);
    return status;
}
