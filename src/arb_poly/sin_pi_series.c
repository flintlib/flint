/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_sin_pi_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        arb_sin_pi(g, h, prec);
        _arb_vec_zero(g + 1, n - 1);
    }
    else if (n == 2)
    {
        arb_t t;
        arb_init(t);
        arb_sin_cos_pi(g, t, h, prec);
        arb_mul(g + 1, h + 1, t, prec);  /* safe since hlen >= 2 */
        arb_const_pi(t, prec);
        arb_mul(g + 1, g + 1, t, prec);
        arb_clear(t);
    }
    else
    {
        arb_ptr t = _arb_vec_init(n);
        _arb_poly_sin_cos_pi_series(g, t, h, hlen, n, prec);
        _arb_vec_clear(t, n);
    }
}

void
arb_poly_sin_pi_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (hlen == 0 || n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    if (hlen == 1)
        n = 1;

    arb_poly_fit_length(g, n);
    _arb_poly_sin_pi_series(g->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}

