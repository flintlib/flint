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

/*
    Third-order Newton iteration: given g = h^(-1/2) + O(x^m), let
    e = 1 - h g^2 = O(x^m). Then

        h^(-1/2) = g (1 - e)^(-1/2) = g (1 + e/2 + 3e^2/8) + O(x^(3m)).

    Compared to the second-order iteration, each step requires one extra
    short squaring but triples the precision, and in practice this is
    faster since the square g^2 and the middle product (h g^2) are
    shared among more new coefficients. When n <= 2m the step reduces
    to the ordinary second-order update.
*/
int
_gr_poly_rsqrt_series_newton(gr_ptr g,
    gr_srcptr h, slong hlen, slong len, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong a[FLINT_BITS];
    slong i, m, n, k;

    hlen = FLINT_MIN(hlen, len);

    if (len == 0)
        return GR_SUCCESS;

    if (hlen == 1)
    {
        status = _gr_poly_rsqrt_series_basecase(g, h, 1, 1, ctx);
        if (len > 1)
            status |= _gr_vec_zero(GR_ENTRY(g, 1, sz), len - 1, ctx);
        return status;
    }

    cutoff = FLINT_MAX(cutoff, 2);

    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 2) / 3);

    status |= _gr_poly_rsqrt_series_basecase(g, h, FLINT_MIN(hlen, n), n, ctx);

    if (status != GR_SUCCESS)
        return status;

    if (len > n)
    {
        gr_ptr t, u, v;
        slong tlen, hnlen, plen, ualloc, valloc;

        /* Hack: until all rings have a good mulmid, implement both with and without. */
        int have_mulmid = (ctx->methods[GR_METHOD_POLY_MULMID] != (gr_funcptr) _gr_poly_mulmid_generic);

        ualloc = have_mulmid ? len - (len + 2) / 3 : len;
        valloc = (len + 2) / 3;

        GR_TMP_INIT_VEC(t, len + ualloc + valloc, ctx);
        u = GR_ENTRY(t, len, sz);
        v = GR_ENTRY(u, ualloc, sz);

        for (i--; i >= 0; i--)
        {
            gr_ptr w;

            m = n;
            n = a[i];
            k = n - 2 * m;

            tlen = FLINT_MIN(2 * m - 1, n);
            status |= _gr_poly_mullow(t, g, m, g, m, tlen, ctx);

            /* w = (h g^2)[m, n) = -e / x^m */
            hnlen = FLINT_MIN(hlen, n);
            plen = FLINT_MIN(n, tlen + hnlen - 1);

            if (have_mulmid)
            {
                status |= _gr_poly_mulmid(u, t, tlen, h, hnlen, m, plen, ctx);
                w = u;
            }
            else
            {
                status |= _gr_poly_mullow(u, t, tlen, h, hnlen, plen, ctx);
                w = GR_ENTRY(u, m, sz);
            }

            if (plen < n)
                status |= _gr_vec_zero(GR_ENTRY(w, plen - m, sz), n - plen, ctx);

            if (k > 0)
            {
                /* e/2 + 3e^2/8 = -(4w - 3 x^m w^2) / 8 */
                status |= _gr_poly_mullow(v, w, k, w, k, k, ctx);
                status |= _gr_vec_mul_scalar_si(v, v, k, 3, ctx);
                status |= _gr_vec_mul_scalar_2exp_si(w, w, n - m, 2, ctx);
                status |= _gr_vec_sub(GR_ENTRY(w, m, sz), GR_ENTRY(w, m, sz), v, k, ctx);
            }

            status |= _gr_poly_mullow(GR_ENTRY(g, m, sz), g, FLINT_MIN(m, n - m), w, n - m, n - m, ctx);
            status |= _gr_vec_mul_scalar_2exp_si(GR_ENTRY(g, m, sz), GR_ENTRY(g, m, sz), n - m, (k > 0) ? -3 : -1, ctx);
            status |= _gr_vec_neg(GR_ENTRY(g, m, sz), GR_ENTRY(g, m, sz), n - m, ctx);
        }

        GR_TMP_CLEAR_VEC(t, len + ualloc + valloc, ctx);
    }

    return status;
}

int
gr_poly_rsqrt_series_newton(gr_poly_t res, const gr_poly_t h, slong len, slong cutoff, gr_ctx_t ctx)
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
        status = gr_poly_rsqrt_series_newton(t, h, len, cutoff, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_rsqrt_series_newton(res->coeffs, h->coeffs, h->length, len, cutoff, ctx);
    _gr_poly_set_length_normalise(res, len, ctx);
    return status;
}
