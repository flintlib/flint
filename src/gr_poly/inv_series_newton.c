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
_gr_poly_inv_series_newton(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, m, n, Qnlen, Wlen, W2len;
    gr_ptr W;
    slong a[FLINT_BITS];

    if (len == 0)
        return GR_SUCCESS;

    if (Qlen == 0)
        return GR_DOMAIN;

    Qlen = FLINT_MIN(Qlen, len);

    if (len < cutoff)
        return _gr_poly_inv_series_basecase(Qinv, Q, Qlen, len, ctx);

    cutoff = FLINT_MAX(cutoff, 2);

    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 1) / 2);

    status |= _gr_poly_inv_series_basecase(Qinv, Q, Qlen, n, ctx);

    if (status != GR_SUCCESS)
        return status;

    GR_TMP_INIT_VEC(W, len, ctx);

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];

        Qnlen = FLINT_MIN(Qlen, n);
        Wlen = FLINT_MIN(Qnlen + m - 1, n);
        W2len = Wlen - m;
        status |= _gr_poly_mullow(W, Q, Qnlen, Qinv, m, Wlen, ctx);
        /* should be mulmid */
        status |= _gr_poly_mullow(GR_ENTRY(Qinv, m, sz), Qinv, m, GR_ENTRY(W, m, sz), W2len, n - m, ctx);
        status |= _gr_vec_neg(GR_ENTRY(Qinv, m, sz), GR_ENTRY(Qinv, m, sz), n - m, ctx);
    }

    GR_TMP_CLEAR_VEC(W, len, ctx);

    return status;
}

int
gr_poly_inv_series_newton(gr_poly_t Qinv, const gr_poly_t Q, slong len, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong Qlen;

    if (len == 0)
        return gr_poly_zero(Qinv, ctx);

    Qlen = Q->length;

    if (Qlen == 0)
        return GR_DOMAIN;

    if (Qlen == 1)
        len = 1;

    if (Qinv == Q)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_inv_series_newton(t, Q, len, cutoff, ctx);
        gr_poly_swap(Qinv, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Qinv, len, ctx);
    status |= _gr_poly_inv_series_newton(Qinv->coeffs, Q->coeffs, Q->length, len, cutoff, ctx);
    _gr_poly_set_length(Qinv, len, ctx);
    _gr_poly_normalise(Qinv, ctx);
    return status;
}
