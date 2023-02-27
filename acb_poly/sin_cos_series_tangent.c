/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_sin_cos_series_tangent(acb_ptr s, acb_ptr c,
        const acb_srcptr h, slong hlen, slong len, slong prec, int times_pi)
{
    acb_ptr t, u, v;
    acb_t s0, c0;
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        if (times_pi)
            acb_sin_cos_pi(s, c, h, prec);
        else
            acb_sin_cos(s, c, h, prec);
        _acb_vec_zero(s + 1, len - 1);
        _acb_vec_zero(c + 1, len - 1);
        return;
    }

    /*
    sin(x) = 2*tan(x/2)/(1+tan(x/2)^2)
    cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2)
    */

    acb_init(s0);
    acb_init(c0);

    t = _acb_vec_init(3 * len);
    u = t + len;
    v = u + len;

    /* sin, cos of h0 */
    if (times_pi)
        acb_sin_cos_pi(s0, c0, h, prec);
    else
        acb_sin_cos(s0, c0, h, prec);

    /* t = tan((h-h0)/2) */
    acb_zero(u);
    _acb_vec_scalar_mul_2exp_si(u + 1, h + 1, hlen - 1, -1);
    if (times_pi)
    {
        acb_const_pi(t, prec);
        _acb_vec_scalar_mul(u + 1, u + 1, hlen - 1, t, prec);
    }

    _acb_poly_tan_series(t, u, hlen, len, prec);

    /* v = 1 + t^2 */
    _acb_poly_mullow(v, t, len, t, len, len, prec);
    acb_add_ui(v, v, 1, prec);

    /* u = 1/(1+t^2) */
    _acb_poly_inv_series(u, v, len, len, prec);

    /* sine */
    _acb_poly_mullow(s, t, len, u, len, len, prec);
    _acb_vec_scalar_mul_2exp_si(s, s, len, 1);

    /* cosine */
    acb_sub_ui(v, v, 2, prec);
    _acb_vec_neg(v, v, len);
    _acb_poly_mullow(c, v, len, u, len, len, prec);

    /* sin(h0 + h1) = cos(h0) sin(h1) + sin(h0) cos(h1)
       cos(h0 + h1) = cos(h0) cos(h1) - sin(h0) sin(h1) */
    if (!acb_is_zero(s0))
    {
        _acb_vec_scalar_mul(t, s, len, c0, prec);
        _acb_vec_scalar_mul(u, c, len, s0, prec);
        _acb_vec_scalar_mul(v, s, len, s0, prec);
        _acb_vec_add(s, t, u, len, prec);
        _acb_vec_scalar_mul(t, c, len, c0, prec);
        _acb_vec_sub(c, t, v, len, prec);
    }

    _acb_vec_clear(t, 3 * len);

    acb_clear(s0);
    acb_clear(c0);
}

void
acb_poly_sin_cos_series_tangent(acb_poly_t s, acb_poly_t c,
                                    const acb_poly_t h, slong n, slong prec, int times_pi)
{
    slong hlen = h->length;

    if (n == 0)
    {
        acb_poly_zero(s);
        acb_poly_zero(c);
        return;
    }

    if (hlen == 0)
    {
        acb_poly_zero(s);
        acb_poly_one(c);
        return;
    }

    acb_poly_fit_length(s, n);
    acb_poly_fit_length(c, n);
    _acb_poly_sin_cos_series_tangent(s->coeffs, c->coeffs,
        h->coeffs, hlen, n, prec, times_pi);
    _acb_poly_set_length(s, n);
    _acb_poly_normalise(s);
    _acb_poly_set_length(c, n);
    _acb_poly_normalise(c);
}

