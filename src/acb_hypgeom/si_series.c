/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_hypgeom_si_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec)
{
    acb_t c;

    acb_init(c);
    acb_hypgeom_si(c, h, prec);

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        _acb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        acb_sinc(g, h, prec);
        acb_mul(g + 1, g, h + 1, prec);
    }
    else
    {
        acb_ptr t, u;

        t = _acb_vec_init(len - 1);
        u = _acb_vec_init(hlen - 1);

        /* Si(h(x)) = integral(h'(x) sinc(h(x))) */
        _acb_poly_sinc_series(t, h, hlen, len - 1, prec);
        _acb_poly_derivative(u, h, hlen, prec);
        _acb_poly_mullow(g, t, len - 1, u, hlen - 1, len - 1, prec);
        _acb_poly_integral(g, g, len, prec);

        _acb_vec_clear(t, len - 1);
        _acb_vec_clear(u, hlen - 1);
    }

    acb_swap(g, c);
    acb_clear(c);
}

void
acb_hypgeom_si_series(acb_poly_t g, const acb_poly_t h, slong len, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || len == 0)
    {
        acb_poly_zero(g);
        return;
    }

    acb_poly_fit_length(g, len);
    _acb_hypgeom_si_series(g->coeffs, h->coeffs, hlen, len, prec);
    _acb_poly_set_length(g, len);
    _acb_poly_normalise(g);
}

