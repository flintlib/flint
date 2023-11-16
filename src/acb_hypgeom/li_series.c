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
_acb_hypgeom_li_series(acb_ptr g, acb_srcptr h, slong hlen, int offset, slong len, slong prec)
{
    acb_t c;

    if (acb_contains_zero(h))
    {
        _acb_vec_indeterminate(g, len);
        return;
    }

    acb_init(c);
    acb_hypgeom_li(c, h, offset, prec);

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        _acb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        acb_log(g, h, prec);
        acb_div(g + 1, h + 1, g, prec);
    }
    else
    {
        acb_ptr t, u;

        t = _acb_vec_init(len);
        u = _acb_vec_init(hlen);

        /* li(h(x)) = integral(h'(x) / log(h(x))) */
        _acb_poly_log_series(t, h, hlen, len - 1, prec);
        _acb_poly_derivative(u, h, hlen, prec);
        _acb_poly_div_series(g, u, hlen - 1, t, len - 1, len - 1, prec);
        _acb_poly_integral(g, g, len, prec);

        _acb_vec_clear(t, len);
        _acb_vec_clear(u, hlen);
    }

    acb_swap(g, c);
    acb_clear(c);
}

void
acb_hypgeom_li_series(acb_poly_t g, const acb_poly_t h, int offset, slong len, slong prec)
{
    slong hlen = h->length;

    if (len == 0)
    {
        acb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        acb_poly_inv_series(g, h, len, prec);
        return;
    }

    acb_poly_fit_length(g, len);
    _acb_hypgeom_li_series(g->coeffs, h->coeffs, hlen, offset, len, prec);
    _acb_poly_set_length(g, len);
    _acb_poly_normalise(g);
}

