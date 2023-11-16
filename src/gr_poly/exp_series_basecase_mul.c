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

/* exp(a + b*x^m) = exp(a) * (1 + b + b^2*x^(2m)/2) */
/* idea: for the Newton basecase, simultaneously compute exp(-a) * (1 + b - b^2*x^(2m)) ? */

int
_gr_poly_exp_series_basecase_mul(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    slong m, v;
    gr_ptr t, u;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    hlen = FLINT_MIN(n, hlen);

    m = (n + 2) / 3;
    v = m * 2;

    if (hlen - m < 1 || n - v < 1 || hlen - v < 1)
        return _gr_poly_exp_series_basecase(f, h, hlen, n, ctx);

    GR_TMP_INIT_VEC(t, 2 * n - m, ctx);
    u = GR_ENTRY(t, n, sz);

    status |= _gr_poly_mullow(t, GR_ENTRY(h, m, sz), hlen - m, GR_ENTRY(h, m, sz), hlen - m, n - v, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(t, t, n - v, -1, ctx);
    status |= _gr_vec_set(u, GR_ENTRY(h, m, sz), v - m, ctx);
    status |= _gr_poly_add(GR_ENTRY(u, v - m, sz), t, n - v, GR_ENTRY(h, v, sz), hlen - v, ctx);
    status |= _gr_poly_exp_series_basecase(f, h, m, n, ctx);
    status |= _gr_poly_mullow(t, f, n, u, n - m, n - m, ctx);
    status |= _gr_poly_add(GR_ENTRY(f, m, sz), GR_ENTRY(f, m, sz), n - m, t, n - m, ctx);

    GR_TMP_CLEAR_VEC(t, 2 * n - m, ctx);

    return status;
}

int
gr_poly_exp_series_basecase_mul(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen = h->length;

    if (n == 0)
        return gr_poly_zero(f, ctx);

    if (hlen == 0)
        return gr_poly_one(f, ctx);

    gr_poly_fit_length(f, n, ctx);
    status |= _gr_poly_exp_series_basecase_mul(f->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(f, n, ctx);
    _gr_poly_normalise(f, ctx);
    return status;
}
