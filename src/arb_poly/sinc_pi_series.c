/*
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_sinc_pi_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        arb_sinc_pi(g, h, prec);
        _arb_vec_zero(g + 1, n - 1);
    }
    else
    {
        arb_t pi;
        arb_ptr t, u;

        arb_init(pi);
        t = _arb_vec_init(n + 1);
        u = _arb_vec_init(hlen);

        arb_const_pi(pi, prec);
        _arb_vec_set(u, h, hlen);

        if (arb_is_zero(h))
        {
            _arb_poly_sin_pi_series(t, u, hlen, n + 1, prec);
            _arb_poly_div_series(g, t + 1, n, u + 1, hlen - 1, n, prec);
        }
        else
        {
            _arb_poly_sin_pi_series(t, u, hlen, n, prec);
            _arb_poly_div_series(g, t, n, u, hlen, n, prec);
        }
        _arb_vec_scalar_div(g, g, n, pi, prec);

        arb_clear(pi);
        _arb_vec_clear(t, n + 1);
        _arb_vec_clear(u, hlen);
    }
}

void
arb_poly_sinc_pi_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_one(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    arb_poly_fit_length(g, n);
    _arb_poly_sinc_pi_series(g->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}

