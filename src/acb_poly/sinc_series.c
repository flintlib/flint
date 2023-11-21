/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_sinc_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        acb_sinc(g, h, prec);
        _acb_vec_zero(g + 1, n - 1);
    }
    else
    {
        acb_ptr t, u;

        t = _acb_vec_init(n + 1);
        u = _acb_vec_init(hlen);

        _acb_vec_set(u, h, hlen);

        if (acb_is_zero(h))
        {
            _acb_poly_sin_series(t, u, hlen, n + 1, prec);
            _acb_poly_div_series(g, t + 1, n, u + 1, hlen - 1, n, prec);
        }
        else
        {
            _acb_poly_sin_series(t, u, hlen, n, prec);
            _acb_poly_div_series(g, t, n, u, hlen, n, prec);
        }

        _acb_vec_clear(t, n + 1);
        _acb_vec_clear(u, hlen);
    }
}

void
acb_poly_sinc_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        acb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        acb_poly_one(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    acb_poly_fit_length(g, n);
    _acb_poly_sinc_series(g->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}

