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

/* c_k x^k -> c_k x^k / (m+k) */
/* todo: optimize */
int
_gr_poly_integral_offset(gr_ptr res, gr_srcptr poly, slong len, slong m, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong k;
    gr_ptr t;
    slong sz = ctx->sizeof_elem;

    if (gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
    {
        GR_TMP_INIT(t, ctx);

        status |= gr_one(t, ctx);
        for (k = len - 1; k >= 0; k--)
        {
            status |= gr_mul(GR_ENTRY(res, k, sz), GR_ENTRY(poly, k, sz), t, ctx);
            status |= gr_mul_ui(t, t, m + k, ctx);
        }

        status |= gr_inv(t, t, ctx);
        /* early return? */

        for (k = 0; k < len; k++)
        {
            status |= gr_mul(GR_ENTRY(res, k, sz), GR_ENTRY(res, k, sz), t, ctx);
            status |= gr_mul_ui(t, t, m + k, ctx);
        }

        GR_TMP_CLEAR(t, ctx);
    }
    else
    {
        for (k = 0; k < len; k++)
            status |= gr_div_ui(GR_ENTRY(res, k, sz), GR_ENTRY(poly, k, sz), m + k, ctx);
    }

    return status;
}

int
_gr_poly_exp_series_newton(gr_ptr f, gr_ptr g,
    gr_srcptr h, slong hlen, slong len, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong a[FLINT_BITS];
    slong n, i, m, l, r;
    gr_ptr t, hprime;
    int inverse;

    cutoff = FLINT_MAX(cutoff, 2);
    hlen = FLINT_MIN(hlen, len);

    if (hlen <= 1 || len <= 1)
    {
        status |= _gr_poly_exp_series_basecase(f, h, hlen, len, ctx);
        if (g != NULL)
            status |= _gr_poly_inv_series(g, f, len, len, ctx);
        return status;
    }

    /* If g is provided, we compute g = exp(-h), and we can use g as
       scratch space. Otherwise, we still need to compute exp(-h) to length
       (n+1)/2 for intermediate use, and we still need n coefficients of
       scratch space. */
    inverse = (g != NULL);
    if (!inverse)
    {
        GR_TMP_INIT_VEC(g, len, ctx);
    }

    hlen = FLINT_MIN(hlen, len);

    GR_TMP_INIT_VEC(t, len, ctx);
    GR_TMP_INIT_VEC(hprime, hlen - 1, ctx);

    status |= _gr_poly_derivative(hprime, h, hlen, ctx);

    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 1) / 2);

    /* f := exp(h) + O(x^n),  g := exp(-h) + O(x^n) */
    status |= _gr_poly_exp_series_basecase_mul(f, h, hlen, n, ctx);
    status |= _gr_poly_inv_series(g, f, n, n, ctx);

    for (i--; i >= 0; i--)
    {
        m = n;             /* previous length */
        n = a[i];          /* new length */

        l = FLINT_MIN(hlen, n) - 1;
        r = FLINT_MIN(l + m - 1, n - 1);
        status |= _gr_poly_mullow(t, hprime, l, f, m, r, ctx);
        status |= _gr_poly_mullow(GR_ENTRY(g, m, sz), g, n - m, GR_ENTRY(t, m - 1, sz), r + 1 - m, n - m, ctx);
        status |= _gr_poly_integral_offset(GR_ENTRY(g, m, sz), GR_ENTRY(g, m, sz), n - m, m, ctx);
        status |= _gr_poly_mullow(GR_ENTRY(f, m, sz), f, n - m, GR_ENTRY(g, m, sz), n - m, n - m, ctx);

        /* g := exp(-h) + O(x^n); not needed if we only want exp(x) */
        if (i != 0 || inverse)
        {
            status |= _gr_poly_mullow(t, f, n, g, m, n, ctx);
            status |= _gr_poly_mullow(GR_ENTRY(g, m, sz), g, m, GR_ENTRY(t, m, sz), n - m, n - m, ctx);
            status |= _gr_vec_neg(GR_ENTRY(g, m, sz), GR_ENTRY(g, m, sz), n - m, ctx);
        }
    }

    GR_TMP_CLEAR_VEC(hprime, hlen - 1, ctx);
    GR_TMP_CLEAR_VEC(t, len, ctx);

    if (!inverse)
    {
        GR_TMP_CLEAR_VEC(g, len, ctx);
    }

    return status;
}

int
gr_poly_exp_series_newton(gr_poly_t f, const gr_poly_t h, slong n, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen = h->length;

    if (n == 0)
        return gr_poly_zero(f, ctx);

    if (hlen == 0)
        return gr_poly_one(f, ctx);

    gr_poly_fit_length(f, n, ctx);
    status |= _gr_poly_exp_series_newton(f->coeffs, NULL, h->coeffs, hlen, n, cutoff, ctx);
    _gr_poly_set_length(f, n, ctx);
    _gr_poly_normalise(f, ctx);
    return status;
}
