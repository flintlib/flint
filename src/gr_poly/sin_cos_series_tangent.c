/*
    Copyright (C) 2013 Fredrik Johansson

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
_gr_poly_sin_cos_series_tangent(gr_ptr s, gr_ptr c,
        gr_srcptr h, slong hlen, slong len, int times_pi, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr t, u, v, s0, c0;

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        if (times_pi)
            status |= gr_sin_cos_pi(s, c, h, ctx);
        else
            status |= gr_sin_cos(s, c, h, ctx);
        status |= _gr_vec_zero(GR_ENTRY(s, 1, sz), len - 1, ctx);
        status |= _gr_vec_zero(GR_ENTRY(c, 1, sz), len - 1, ctx);
        return status;
    }

    /*
    sin(x) = 2*tan(x/2)/(1+tan(x/2)^2)
    cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2)
    */

    GR_TMP_INIT_VEC(t, 3 * len + 2, ctx);
    u = GR_ENTRY(t, len, sz);
    v = GR_ENTRY(u, len, sz);
    s0 = GR_ENTRY(v, len, sz);
    c0 = GR_ENTRY(s0, 1, sz);

    /* sin, cos of h0 */
    if (times_pi)
        status |= gr_sin_cos_pi(s0, c0, h, ctx);
    else
        status |= gr_sin_cos(s0, c0, h, ctx);

    /* t = tan((h-h0)/2) */
    status |= gr_zero(u, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(GR_ENTRY(u, 1, sz), GR_ENTRY(h, 1, sz), hlen - 1, -1, ctx);
    if (times_pi)
    {
        status |= gr_pi(t, ctx);
        status |= _gr_vec_mul_scalar(GR_ENTRY(u, 1, sz), GR_ENTRY(u, 1, sz), hlen - 1, t, ctx);
    }

    status |= _gr_poly_tan_series(t, u, hlen, len, ctx);

    /* v = 1 + t^2 */
    status |= _gr_poly_mullow(v, t, len, t, len, len, ctx);
    status |= gr_add_ui(v, v, 1, ctx);

    /* u = 1/(1+t^2) */
    status |= _gr_poly_inv_series(u, v, len, len, ctx);

    /* sine */
    status |= _gr_poly_mullow(s, t, len, u, len, len, ctx);
    status |= _gr_vec_mul_scalar_2exp_si(s, s, len, 1, ctx);

    /* cosine */
    status |= gr_sub_ui(v, v, 2, ctx);
    status |= _gr_vec_neg(v, v, len, ctx);
    status |= _gr_poly_mullow(c, v, len, u, len, len, ctx);

    /* sin(h0 + h1) = cos(h0) sin(h1) + sin(h0) cos(h1)
       cos(h0 + h1) = cos(h0) cos(h1) - sin(h0) sin(h1) */
    if (gr_is_zero(s0, ctx) != T_TRUE)
    {
        status |= _gr_vec_mul_scalar(t, s, len, c0, ctx);
        status |= _gr_vec_mul_scalar(u, c, len, s0, ctx);
        status |= _gr_vec_mul_scalar(v, s, len, s0, ctx);
        status |= _gr_vec_add(s, t, u, len, ctx);
        status |= _gr_vec_mul_scalar(t, c, len, c0, ctx);
        status |= _gr_vec_sub(c, t, v, len, ctx);
    }

    GR_TMP_CLEAR_VEC(t, 3 * len + 2, ctx);

    return status;
}

int
gr_poly_sin_cos_series_tangent(gr_poly_t s, gr_poly_t c,
                                    const gr_poly_t h, slong n, int times_pi, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen = h->length;

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

    gr_poly_fit_length(s, n, ctx);
    gr_poly_fit_length(c, n, ctx);
    status |= _gr_poly_sin_cos_series_tangent(s->coeffs, c->coeffs,
        h->coeffs, hlen, n, times_pi, ctx);
    _gr_poly_set_length(s, n, ctx);
    _gr_poly_normalise(s, ctx);
    _gr_poly_set_length(c, n, ctx);
    _gr_poly_normalise(c, ctx);

    return status;
}
