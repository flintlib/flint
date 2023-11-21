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

int
_gr_poly_sqrt_series_newton(gr_ptr g,
    gr_srcptr h, slong hlen, slong len, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong a[FLINT_BITS];
    slong i, m, n, alloc;
    gr_ptr t, u, v;
    slong tlen, ulen;

    hlen = FLINT_MIN(hlen, len);

    if (len == 0)
        return GR_SUCCESS;

    if (len < cutoff)
        return _gr_poly_sqrt_series_basecase(g, h, hlen, len, ctx);

    cutoff = FLINT_MAX(cutoff, 2);
    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 1) / 2);

    status |= _gr_poly_rsqrt_series_basecase(g, h, FLINT_MIN(hlen, n), n, ctx);

    if (status != GR_SUCCESS)
        return status;

    alloc = 2 * len + (len + 1) / 2;

    GR_TMP_INIT_VEC(t, alloc, ctx);
    u = GR_ENTRY(t, len, sz);
    v = GR_ENTRY(u, len, sz);

    for (i--; i >= 1; i--)
    {
        m = n;
        n = a[i];

        tlen = FLINT_MIN(2 * m - 1, n);
        ulen = FLINT_MIN(n, m + tlen - 1);

        status |= _gr_poly_mullow(t, g, m, g, m, tlen, ctx);
        status |= _gr_poly_mullow(u, g, m, t, tlen, ulen, ctx);
        status |= _gr_poly_mullow(t, u, ulen, h, FLINT_MIN(hlen, n), n, ctx);   /* should be mulmid */
        status |= _gr_vec_mul_scalar_2exp_si(GR_ENTRY(g, m, sz), GR_ENTRY(t, m, sz), n - m, -1, ctx);
        status |= _gr_vec_neg(GR_ENTRY(g, m, sz), GR_ENTRY(g, m, sz), n - m, ctx);
    }

    m = (len + 1) / 2;
    n = len;

    /* Karp-Markstein */
    /* todo: cleanup; improve allocations? */
    tlen = FLINT_MIN(2 * m - 1, n);
    status |= _gr_poly_mullow(v, g, m, h, hlen, m, ctx);
    status |= _gr_poly_mullow(t, v, m, v, m, tlen, ctx);
    status |= _gr_poly_sub(GR_ENTRY(u, m, sz), GR_ENTRY(h, m, sz), FLINT_MAX(0, FLINT_MIN(hlen - m, n - m)), GR_ENTRY(t, m, sz), FLINT_MAX(0, FLINT_MIN(tlen - m, n - m)), ctx);
    status |= _gr_poly_mullow(t, g, m, GR_ENTRY(u, m, sz), n - m, n - m, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(GR_ENTRY(g, m, sz), t, n - m, -1, ctx);
    _gr_vec_swap(g, v, m, ctx);

    GR_TMP_CLEAR_VEC(t, alloc, ctx);

    return status;
}

int
gr_poly_sqrt_series_newton(gr_poly_t res, const gr_poly_t h, slong len, slong cutoff, gr_ctx_t ctx)
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
