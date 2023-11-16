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
_arb_hypgeom_si_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec)
{
    arb_t c;

    arb_init(c);
    arb_hypgeom_si(c, h, prec);

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        _arb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        arb_sinc(g, h, prec);
        arb_mul(g + 1, g, h + 1, prec);
    }
    else
    {
        arb_ptr t, u;

        t = _arb_vec_init(len - 1);
        u = _arb_vec_init(hlen - 1);

        /* Si(h(x)) = integral(h'(x) sinc(h(x))) */
        _arb_poly_sinc_series(t, h, hlen, len - 1, prec);
        _arb_poly_derivative(u, h, hlen, prec);
        _arb_poly_mullow(g, t, len - 1, u, hlen - 1, len - 1, prec);
        _arb_poly_integral(g, g, len, prec);

        _arb_vec_clear(t, len - 1);
        _arb_vec_clear(u, hlen - 1);
    }

    arb_swap(g, c);
    arb_clear(c);
}

void
arb_hypgeom_si_series(arb_poly_t g, const arb_poly_t h, slong len, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || len == 0)
    {
        arb_poly_zero(g);
        return;
    }

    arb_poly_fit_length(g, len);
    _arb_hypgeom_si_series(g->coeffs, h->coeffs, hlen, len, prec);
    _arb_poly_set_length(g, len);
    _arb_poly_normalise(g);
}

