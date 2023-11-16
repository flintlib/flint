/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_hypgeom.h"

void
_arb_hypgeom_gamma_upper_series(arb_ptr g, const arb_t s, arb_srcptr h, slong hlen, int regularized, slong n, slong prec)
{
    arb_t c;
    arb_init(c);

    arb_hypgeom_gamma_upper(c, s, h, regularized, prec);

    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        _arb_vec_zero(g + 1, n - 1);
    }
    else
    {
        arb_ptr t, u, v;
        arb_ptr w = NULL;

        t = _arb_vec_init(n);
        u = _arb_vec_init(n);
        v = _arb_vec_init(n);

        if (regularized == 2)
        {
            w = _arb_vec_init(n);
            arb_neg(t, s);
            _arb_poly_pow_arb_series(w, h, hlen, t, n, prec);
        }

        /* Gamma(s, h(x)) = -integral(h'(x) h(x)^(s-1) exp(-h(x)) */
        arb_sub_ui(u, s, 1, prec);
        _arb_poly_pow_arb_series(t, h, hlen, u, n, prec);
        _arb_poly_derivative(u, h, hlen, prec);
        _arb_poly_mullow(v, t, n, u, hlen - 1, n, prec);
        _arb_vec_neg(t, h, hlen);
        _arb_poly_exp_series(t, t, hlen, n, prec);
        _arb_poly_mullow(g, v, n, t, n, n, prec);
        _arb_poly_integral(g, g, n, prec);
        _arb_vec_neg(g, g, n);

        if (regularized == 1)
        {
            arb_rgamma(t, s, prec);
            _arb_vec_scalar_mul(g, g, n, t, prec);
        }
        else if (regularized == 2)
        {
            _arb_vec_set(u, g, n);
            _arb_poly_mullow(g, u, n, w, n, n, prec);
            _arb_vec_clear(w, n);
        }

        _arb_vec_clear(t, n);
        _arb_vec_clear(u, n);
        _arb_vec_clear(v, n);
    }

    arb_swap(g, c);
    arb_clear(c);
}

void
arb_hypgeom_gamma_upper_series(arb_poly_t g, const arb_t s,
        const arb_poly_t h, int regularized, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    arb_poly_fit_length(g, n);

    if (hlen == 0)
    {
        arb_t t;
        arb_init(t);
        _arb_hypgeom_gamma_upper_series(g->coeffs, s, t, 1, regularized, n, prec);
        arb_clear(t);
    }
    else
    {
        _arb_hypgeom_gamma_upper_series(
                g->coeffs, s, h->coeffs, hlen, regularized, n, prec);
    }

    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}
