/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_sinh_cosh_series_exponential(acb_ptr s, acb_ptr c,
    const acb_srcptr h, slong hlen, slong len, slong prec)
{
    acb_ptr t, u, v;
    acb_t s0, c0;
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_sinh_cosh(s, c, h, prec);
        _acb_vec_zero(s + 1, len - 1);
        _acb_vec_zero(c + 1, len - 1);
        return;
    }

    acb_init(s0);
    acb_init(c0);

    t = _acb_vec_init(3 * len);
    u = t + len;
    v = u + len;

    acb_sinh_cosh(s0, c0, h, prec);

    _acb_vec_set(t + 1, h + 1, hlen - 1);
    _acb_poly_exp_series(t, t, len, len, prec);

    /* todo: part of the inverse could be avoided since exp computes
       it internally to half the length */
    _acb_poly_inv_series(u, t, len, len, prec);

    /* hyperbolic sine */
    _acb_vec_sub(s, t, u, len, prec);
    _acb_vec_scalar_mul_2exp_si(s, s, len, -1);

    /* hyperbolic cosine */
    _acb_vec_add(c, t, u, len, prec);
    _acb_vec_scalar_mul_2exp_si(c, c, len, -1);

    /* sinh(h0 + h1) = cosh(h0) sinh(h1) + sinh(h0) cosh(h1)
       cosh(h0 + h1) = cosh(h0) cosh(h1) + sinh(h0) sinh(h1) */
    if (!acb_is_zero(s0))
    {
        _acb_vec_scalar_mul(t, s, len, c0, prec);
        _acb_vec_scalar_mul(u, c, len, s0, prec);
        _acb_vec_scalar_mul(v, s, len, s0, prec);
        _acb_vec_add(s, t, u, len, prec);
        _acb_vec_scalar_mul(t, c, len, c0, prec);
        _acb_vec_add(c, t, v, len, prec);
    }

    _acb_vec_clear(t, 3 * len);

    acb_clear(s0);
    acb_clear(c0);
}

void
acb_poly_sinh_cosh_series_exponential(acb_poly_t s, acb_poly_t c,
                                     const acb_poly_t h, slong n, slong prec)
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
    _acb_poly_sinh_cosh_series_exponential(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(s, n);
    _acb_poly_normalise(s);
    _acb_poly_set_length(c, n);
    _acb_poly_normalise(c);
}

