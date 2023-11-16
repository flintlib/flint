/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_hypgeom.h"

void
_arb_hypgeom_ei_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec)
{
    arb_t c;

    if (arb_contains_zero(h))
    {
        _arb_vec_indeterminate(g, len);
        return;
    }

    arb_init(c);
    arb_hypgeom_ei(c, h, prec);

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        _arb_vec_zero(g + 1, len - 1);
    }
    else
    {
        arb_ptr t, u, v;

        t = _arb_vec_init(len);
        u = _arb_vec_init(len);
        v = _arb_vec_init(len);

        /* Ei(h(x)) = integral([h'(x) / h(x)] exp(h(x)) */
        _arb_poly_exp_series(t, h, hlen, len - 1, prec);
        _arb_poly_derivative(u, h, hlen, prec);
        _arb_poly_mullow(v, t, len - 1, u, len - 1, len - 1, prec);
        _arb_poly_div_series(u, v, len - 1, h, hlen, len - 1, prec);
        _arb_poly_integral(g, u, len, prec);

        _arb_vec_clear(t, len);
        _arb_vec_clear(u, len);
        _arb_vec_clear(v, len);
    }

    arb_swap(g, c);
    arb_clear(c);
}

void
arb_hypgeom_ei_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec)
{
    slong hlen = h->length;

    if (len == 0)
    {
        arb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_inv_series(g, h, len, prec);
        return;
    }

    arb_poly_fit_length(g, len);
    _arb_hypgeom_ei_series(g->coeffs, h->coeffs, hlen, len, prec);
    _arb_poly_set_length(g, len);
    _arb_poly_normalise(g);
}

