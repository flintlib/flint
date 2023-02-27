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
_arb_poly_atan_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
{
    arb_t c;
    arb_init(c);

    arb_atan(c, h, prec);

    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        _arb_vec_zero(g + 1, n - 1);
    }
    else
    {
        arb_ptr t, u;
        slong ulen;

        t = _arb_vec_init(n);
        u = _arb_vec_init(n);

        /* atan(h(x)) = integral(h'(x)/(1+h(x)^2)) */
        ulen = FLINT_MIN(n, 2 * hlen - 1);
        _arb_poly_mullow(u, h, hlen, h, hlen, ulen, prec);
        arb_add_ui(u, u, 1, prec);

        _arb_poly_derivative(t, h, hlen, prec);
        _arb_poly_div_series(g, t, hlen - 1, u, ulen, n, prec);
        _arb_poly_integral(g, g, n, prec);

        _arb_vec_clear(t, n);
        _arb_vec_clear(u, n);
    }

    arb_swap(g, c);
    arb_clear(c);
}

void
arb_poly_atan_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    arb_poly_fit_length(g, n);
    _arb_poly_atan_series(g->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}

