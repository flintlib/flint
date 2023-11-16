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
_arb_poly_sqrt_series(arb_ptr g,
    arb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    while (hlen > 0 && arb_is_zero(h + hlen - 1))
        hlen--;

    if (hlen <= 1)
    {
        arb_sqrt(g, h, prec);
        _arb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        arb_sqrt(g, h, prec);
        arb_div(g + 1, h + 1, h, prec);
        arb_mul(g + 1, g + 1, g, prec);
        arb_mul_2exp_si(g + 1, g + 1, -1);
    }
    else if (_arb_vec_is_zero(h + 1, hlen - 2))
    {
        arb_t t;
        arb_init(t);
        arf_set_si_2exp_si(arb_midref(t), 1, -1);
        _arb_poly_binomial_pow_arb_series(g, h, hlen, t, len, prec);
        arb_clear(t);
    }
    else
    {
        arb_ptr t;
        t = _arb_vec_init(len);
        _arb_poly_rsqrt_series(t, h, hlen, len, prec);
        _arb_poly_mullow(g, t, len, h, hlen, len, prec);
        _arb_vec_clear(t, len);
    }
}

void
arb_poly_sqrt_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    if (g == h)
    {
        arb_poly_t t;
        arb_poly_init(t);
        arb_poly_sqrt_series(t, h, n, prec);
        arb_poly_swap(g, t);
        arb_poly_clear(t);
        return;
    }

    arb_poly_fit_length(g, n);
    if (h->length == 0)
        _arb_vec_indeterminate(g->coeffs, n);
    else
        _arb_poly_sqrt_series(g->coeffs, h->coeffs, h->length, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}

