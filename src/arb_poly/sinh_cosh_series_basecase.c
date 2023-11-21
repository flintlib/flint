/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_sinh_cosh_series_basecase(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen,
        slong n, slong prec)
{
    slong k, alen = FLINT_MIN(n, hlen);
    arb_ptr a;
    arb_t t, u;

    arb_sinh_cosh(s, c, h, prec);

    if (hlen == 1)
    {
        _arb_vec_zero(s + 1, n - 1);
        _arb_vec_zero(c + 1, n - 1);
        return;
    }

    arb_init(t);
    arb_init(u);
    a = _arb_vec_init(alen);

    for (k = 1; k < alen; k++)
        arb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        arb_dot(t, NULL, 0, a + 1, 1, s + k - 1, -1, FLINT_MIN(k, hlen - 1), prec);
        arb_dot(u, NULL, 0, a + 1, 1, c + k - 1, -1, FLINT_MIN(k, hlen - 1), prec);
        arb_div_ui(c + k, t, k, prec);
        arb_div_ui(s + k, u, k, prec);
    }

    arb_clear(t);
    arb_clear(u);
    _arb_vec_clear(a, alen);
}

void
arb_poly_sinh_cosh_series_basecase(arb_poly_t s, arb_poly_t c,
        const arb_poly_t h, slong n, slong prec)
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
    _arb_poly_sinh_cosh_series_basecase(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(s, n);
    _arb_poly_normalise(s);
    _arb_poly_set_length(c, n);
    _arb_poly_normalise(c);
}

