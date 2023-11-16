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
_acb_poly_sinh_cosh_series_basecase(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen,
        slong n, slong prec)
{
    slong k, alen = FLINT_MIN(n, hlen);
    acb_ptr a;
    acb_t t, u;

    acb_sinh_cosh(s, c, h, prec);

    if (hlen == 1)
    {
        _acb_vec_zero(s + 1, n - 1);
        _acb_vec_zero(c + 1, n - 1);
        return;
    }

    acb_init(t);
    acb_init(u);
    a = _acb_vec_init(alen);

    for (k = 1; k < alen; k++)
        acb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        acb_dot(t, NULL, 0, a + 1, 1, s + k - 1, -1, FLINT_MIN(k, hlen - 1), prec);
        acb_dot(u, NULL, 0, a + 1, 1, c + k - 1, -1, FLINT_MIN(k, hlen - 1), prec);
        acb_div_ui(c + k, t, k, prec);
        acb_div_ui(s + k, u, k, prec);
    }

    acb_clear(t);
    acb_clear(u);
    _acb_vec_clear(a, alen);
}

void
acb_poly_sinh_cosh_series_basecase(acb_poly_t s, acb_poly_t c,
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
    _acb_poly_sinh_cosh_series_basecase(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(s, n);
    _acb_poly_normalise(s);
    _acb_poly_set_length(c, n);
    _acb_poly_normalise(c);
}
