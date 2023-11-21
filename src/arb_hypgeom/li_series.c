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
_arb_hypgeom_li_series(arb_ptr g, arb_srcptr h, slong hlen, int offset, slong len, slong prec)
{
    arb_t c;

    if (!arb_is_positive(h) || arb_contains_si(h, 1))
    {
        _arb_vec_indeterminate(g, len);
        return;
    }

    arb_init(c);
    arb_hypgeom_li(c, h, offset, prec);

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        _arb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        arb_log(g, h, prec);
        arb_div(g + 1, h + 1, g, prec);
    }
    else
    {
        arb_ptr t, u;

        t = _arb_vec_init(len);
        u = _arb_vec_init(hlen);

        /* li(h(x)) = integral(h'(x) / log(h(x))) */
        _arb_poly_log_series(t, h, hlen, len - 1, prec);
        _arb_poly_derivative(u, h, hlen, prec);
        _arb_poly_div_series(g, u, hlen - 1, t, len - 1, len - 1, prec);
        _arb_poly_integral(g, g, len, prec);

        _arb_vec_clear(t, len);
        _arb_vec_clear(u, hlen);
    }

    arb_swap(g, c);
    arb_clear(c);
}

void
arb_hypgeom_li_series(arb_poly_t g, const arb_poly_t h, int offset, slong len, slong prec)
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
    _arb_hypgeom_li_series(g->coeffs, h->coeffs, hlen, offset, len, prec);
    _arb_poly_set_length(g, len);
    _arb_poly_normalise(g);
}

