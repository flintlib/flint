/*
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_sinc_pi_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        acb_sinc_pi(g, h, prec);
        _acb_vec_zero(g + 1, n - 1);
    }
    else
    {
        acb_t pi;
        acb_ptr t, u;

        acb_init(pi);
        t = _acb_vec_init(n + 1);
        u = _acb_vec_init(hlen);

        acb_const_pi(pi, prec);
        _acb_vec_set(u, h, hlen);

        if (acb_is_zero(h))
        {
            _acb_poly_sin_pi_series(t, u, hlen, n + 1, prec);
            _acb_poly_div_series(g, t + 1, n, u + 1, hlen - 1, n, prec);
            _acb_vec_scalar_div(g, g, n, pi, prec);
        }
        else if (acb_contains_zero(h))
        {
            _acb_vec_scalar_mul(u, h, hlen, pi, prec);
            _acb_poly_sinc_series(g, u, hlen, n, prec);
        }
        else
        {
            _acb_poly_sin_pi_series(t, u, hlen, n, prec);
            _acb_poly_div_series(g, t, n, u, hlen, n, prec);
            _acb_vec_scalar_div(g, g, n, pi, prec);
        }

        acb_clear(pi);
        _acb_vec_clear(t, n + 1);
        _acb_vec_clear(u, hlen);
    }
}

void
acb_poly_sinc_pi_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
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
    _acb_poly_sinc_pi_series(g->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}
