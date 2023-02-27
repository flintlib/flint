/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_sin_cos_series_tangent(arb_ptr s, arb_ptr c,
        arb_srcptr h, slong hlen, slong len, slong prec, int times_pi)
{
    arb_ptr t, u, v;
    arb_t s0, c0;
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        if (times_pi)
            arb_sin_cos_pi(s, c, h, prec);
        else
            arb_sin_cos(s, c, h, prec);
        _arb_vec_zero(s + 1, len - 1);
        _arb_vec_zero(c + 1, len - 1);
        return;
    }

    /*
    sin(x) = 2*tan(x/2)/(1+tan(x/2)^2)
    cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2)
    */

    arb_init(s0);
    arb_init(c0);

    t = _arb_vec_init(3 * len);
    u = t + len;
    v = u + len;

    /* sin, cos of h0 */
    if (times_pi)
        arb_sin_cos_pi(s0, c0, h, prec);
    else
        arb_sin_cos(s0, c0, h, prec);

    /* t = tan((h-h0)/2) */
    arb_zero(u);
    _arb_vec_scalar_mul_2exp_si(u + 1, h + 1, hlen - 1, -1);
    if (times_pi)
    {
        arb_const_pi(t, prec);
        _arb_vec_scalar_mul(u + 1, u + 1, hlen - 1, t, prec);
    }

    _arb_poly_tan_series(t, u, hlen, len, prec);

    /* v = 1 + t^2 */
    _arb_poly_mullow(v, t, len, t, len, len, prec);
    arb_add_ui(v, v, 1, prec);

    /* u = 1/(1+t^2) */
    _arb_poly_inv_series(u, v, len, len, prec);

    /* sine */
    _arb_poly_mullow(s, t, len, u, len, len, prec);
    _arb_vec_scalar_mul_2exp_si(s, s, len, 1);

    /* cosine */
    arb_sub_ui(v, v, 2, prec);
    _arb_vec_neg(v, v, len);
    _arb_poly_mullow(c, v, len, u, len, len, prec);

    /* sin(h0 + h1) = cos(h0) sin(h1) + sin(h0) cos(h1)
       cos(h0 + h1) = cos(h0) cos(h1) - sin(h0) sin(h1) */
    if (!arb_is_zero(s0))
    {
        _arb_vec_scalar_mul(t, s, len, c0, prec);
        _arb_vec_scalar_mul(u, c, len, s0, prec);
        _arb_vec_scalar_mul(v, s, len, s0, prec);
        _arb_vec_add(s, t, u, len, prec);
        _arb_vec_scalar_mul(t, c, len, c0, prec);
        _arb_vec_sub(c, t, v, len, prec);
    }

    _arb_vec_clear(t, 3 * len);

    arb_clear(s0);
    arb_clear(c0);
}

void
arb_poly_sin_cos_series_tangent(arb_poly_t s, arb_poly_t c,
                                    const arb_poly_t h, slong n, slong prec, int times_pi)
{
    slong hlen = h->length;

    if (n == 0)
    {
        arb_poly_zero(s);
        arb_poly_zero(c);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_zero(s);
        arb_poly_one(c);
        return;
    }

    arb_poly_fit_length(s, n);
    arb_poly_fit_length(c, n);
    _arb_poly_sin_cos_series_tangent(s->coeffs, c->coeffs,
        h->coeffs, hlen, n, prec, times_pi);
    _arb_poly_set_length(s, n);
    _arb_poly_normalise(s);
    _arb_poly_set_length(c, n);
    _arb_poly_normalise(c);
}

