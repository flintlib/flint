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

/*
Assumes:
    f[0] = exp(h[0])
    hprime[k] = (k+1) * h[k+1]
*/
int
_gr_poly_exp_series_basecase_rec_precomp1(gr_ptr f, gr_srcptr hprime, slong hlen, slong n, gr_ctx_t ctx)
{
    slong k, l;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    status |= gr_mul(GR_ENTRY(f, 1, sz), f, hprime, ctx);
    if (n == 2)
        return status;

    for (k = 2; k < n && status == GR_SUCCESS; k++)
    {
        l = FLINT_MIN(k, hlen - 1);
        status |= _gr_vec_dot_rev(GR_ENTRY(f, k, sz), NULL, 0, hprime, GR_ENTRY(f, k - l, sz), l, ctx);
        status |= gr_div_ui(GR_ENTRY(f, k, sz), GR_ENTRY(f, k, sz), k, ctx);
    }

    return status;
}

/*
Assumes:
    f[0] = exp(h[0])
    f[k-1] = 1/k, k >= 2
    hprime[k] = (k+1) * h[k+1]
*/
int
_gr_poly_exp_series_basecase_rec_precomp2(gr_ptr f, gr_srcptr hprime, slong hlen, slong n, gr_ctx_t ctx)
{
    slong k, l;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr t;

    status |= gr_mul(GR_ENTRY(f, 1, sz), f, hprime, ctx);
    if (n == 2)
        return status;

    GR_TMP_INIT(t, ctx);

    for (k = 2; k < n && status == GR_SUCCESS; k++)
    {
        l = FLINT_MIN(k, hlen - 1);
        status |= _gr_vec_dot_rev(t, NULL, 0, hprime, GR_ENTRY(f, k - l, sz), l, ctx);
        status |= gr_mul(GR_ENTRY(f, k, sz), GR_ENTRY(f, k, sz), t, ctx);
    }

    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
_gr_poly_exp_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    hlen = FLINT_MIN(hlen, n);

    status = gr_exp(f, h, ctx);
    if (status != GR_SUCCESS)
        return status;

    if (hlen == 1)
        return _gr_vec_zero(GR_ENTRY(f, 1, sz), n - 1, ctx);

    if (n == 1)
        return status;

    if (n == 2)
    {
        status |= gr_mul(GR_ENTRY(f, 1, sz), f, GR_ENTRY(h, 1, sz), ctx);
        return status;
    }

    if (_gr_vec_is_zero(GR_ENTRY(h, 1, sz), hlen - 2, ctx) == T_TRUE) /* h = a + bx^d */
    {
        slong i, j, d = hlen - 1;
        gr_ptr t;

        GR_TMP_INIT(t, ctx);

        status |= gr_set(t, GR_ENTRY(h, d, sz), ctx);

        for (i = 1, j = d; j < n && status == GR_SUCCESS; j += d, i++)
        {
            status |= gr_mul(GR_ENTRY(f, j, sz), GR_ENTRY(f, j - d, sz), t, ctx);
            status |= gr_div_ui(GR_ENTRY(f, j, sz), GR_ENTRY(f, j, sz), i, ctx);
            status |= _gr_vec_zero(GR_ENTRY(f, j - d + 1, sz), hlen - 2, ctx);
        }

        status |= _gr_vec_zero(GR_ENTRY(f, j - d + 1, sz), n - (j - d + 1), ctx);

        GR_TMP_CLEAR(t, ctx);
        return status;
    }

    GR_TMP_INIT_VEC(t, hlen - 1, ctx);
    status |= _gr_poly_derivative(t, h, hlen, ctx);

    if (n < 450 || gr_ctx_is_finite_characteristic(ctx) != T_TRUE)
    {
        status |= _gr_poly_exp_series_basecase_rec_precomp1(f, t, hlen, n, ctx);
    }
    else
    {
        status |= _gr_vec_reciprocals(GR_ENTRY(f, 1, sz), n - 1, ctx);
        status |= _gr_poly_exp_series_basecase_rec_precomp2(f, t, hlen, n, ctx);
    }

    GR_TMP_CLEAR_VEC(t, hlen - 1, ctx);

    return status;
}

int
gr_poly_exp_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen = h->length;

    if (n == 0)
        return gr_poly_zero(f, ctx);

    if (hlen == 0)
        return gr_poly_one(f, ctx);

    gr_poly_fit_length(f, n, ctx);
    status |= _gr_poly_exp_series_basecase(f->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(f, n, ctx);
    _gr_poly_normalise(f, ctx);
    return status;
}
