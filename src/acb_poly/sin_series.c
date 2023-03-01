/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_sin_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        acb_sin(g, h, prec);
        _acb_vec_zero(g + 1, n - 1);
    }
    else if (n == 2)
    {
        acb_t t;
        acb_init(t);
        acb_sin_cos(g, t, h, prec);
        acb_mul(g + 1, h + 1, t, prec);  /* safe since hlen >= 2 */
        acb_clear(t);
    }
    else
    {
        acb_ptr t = _acb_vec_init(n);
        _acb_poly_sin_cos_series(g, t, h, hlen, n, prec);
        _acb_vec_clear(t, n);
    }
}

void
acb_poly_sin_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || n == 0)
    {
        acb_poly_zero(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    acb_poly_fit_length(g, n);
    _acb_poly_sin_series(g->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}

