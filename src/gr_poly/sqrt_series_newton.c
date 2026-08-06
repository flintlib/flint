/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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
    slong tlen, r1, r2, rlen;

    hlen = FLINT_MIN(hlen, len);

    if (len == 0)
        return GR_SUCCESS;

    if (len < cutoff || len == 1)
        return _gr_poly_sqrt_series_basecase(g, h, hlen, len, ctx);

    if (hlen == 1)
    {
        status = _gr_poly_sqrt_series_basecase(g, h, 1, 1, ctx);
        status |= _gr_vec_zero(GR_ENTRY(g, 1, sz), len - 1, ctx);
        return status;
    }

    cutoff = FLINT_MAX(cutoff, 2);
    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 1) / 2);

    status |= _gr_poly_rsqrt_series_basecase(g, h, FLINT_MIN(hlen, n), n, ctx);

    if (status != GR_SUCCESS)
        return status;

    /* Hack: until all rings have a good mulmid, implement both with and without. */
    int have_mulmid = (ctx->methods[GR_METHOD_POLY_MULMID] != (gr_funcptr) _gr_poly_mulmid_generic);

    alloc = (have_mulmid ? (len + 1) / 2 : len) + 2 * ((len + 1) / 2);

    GR_TMP_INIT_VEC(t, alloc, ctx);
    u = GR_ENTRY(t, have_mulmid ? (len + 1) / 2 : len, sz);
    v = GR_ENTRY(u, (len + 1) / 2, sz);

    for (i--; i >= 1; i--)
    {
        m = n;
        n = a[i];

        tlen = FLINT_MIN(2 * m - 1, n);

        status |= _gr_poly_mullow(t, g, m, g, m, tlen, ctx);

        if (have_mulmid)
        {
            status |= _gr_poly_mulmid(u, t, tlen, h, FLINT_MIN(hlen, n), m, n, ctx);
            status |= _gr_poly_mullow(GR_ENTRY(g, m, sz), g, n - m, u, n - m, n - m, ctx);
        }
        else
        {
            status |= _gr_poly_mullow(u, t, tlen, h, FLINT_MIN(hlen, n), n, ctx);
            status |= _gr_poly_mullow(GR_ENTRY(g, m, sz), g, n - m, GR_ENTRY(u, m, sz), n - m, n - m, ctx);
        }

        status |= _gr_vec_mul_scalar_2exp_si(GR_ENTRY(g, m, sz), GR_ENTRY(g, m, sz), n - m, -1, ctx);
        status |= _gr_vec_neg(GR_ENTRY(g, m, sz), GR_ENTRY(g, m, sz), n - m, ctx);
    }

    m = (len + 1) / 2;
    n = len;

    /* Karp-Markstein */
    tlen = FLINT_MIN(2 * m - 1, n);

    status |= _gr_poly_mullow(v, g, m, h, FLINT_MIN(hlen, m), m, ctx);

    r1 = FLINT_MAX(0, FLINT_MIN(hlen - m, n - m));
    r2 = FLINT_MAX(0, FLINT_MIN(tlen - m, n - m));
    rlen = FLINT_MAX(r1, r2);

    if (have_mulmid)
    {
        if (m < tlen)
            status |= _gr_poly_mulmid(t, v, m, v, m, m, tlen, ctx);
        status |= _gr_poly_sub(u, GR_ENTRY(h, m, sz), r1, t, r2, ctx);
    }
    else
    {
        status |= _gr_poly_mullow(t, v, m, v, m, tlen, ctx);
        status |= _gr_poly_sub(u, GR_ENTRY(h, m, sz), r1, GR_ENTRY(t, m, sz), r2, ctx);
    }

    status |= _gr_poly_mullow(GR_ENTRY(g, m, sz), g, n - m, u, rlen, n - m, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(GR_ENTRY(g, m, sz), GR_ENTRY(g, m, sz), n - m, -1, ctx);
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
