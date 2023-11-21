/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_hypgeom_gamma_lower_series(acb_ptr g, const acb_t s, acb_srcptr h, slong hlen, int regularized, slong n, slong prec)
{
    acb_t c;

    hlen = FLINT_MIN(hlen, n);

    if (regularized == 2 && acb_is_int(s) && arb_is_nonpositive(acb_realref(s)))
    {
        acb_t ns;
        acb_init(ns);
        acb_neg(ns, s);
        if (g == h)
        {
            acb_ptr t = _acb_vec_init(hlen);
            _acb_vec_set(t, h, hlen);
            _acb_poly_pow_acb_series(g, t, hlen, ns, n, prec);
            _acb_vec_clear(t, hlen);
        }
        else
        {
            _acb_poly_pow_acb_series(g, h, hlen, ns, n, prec);
        }
        acb_clear(ns);
        return;
    }

    acb_init(c);
    acb_hypgeom_gamma_lower(c, s, h, regularized, prec);

    if (hlen == 1)
    {
        _acb_vec_zero(g + 1, n - 1);
    }
    else
    {
        acb_ptr t, u, v;
        acb_ptr w = NULL;

        t = _acb_vec_init(n);
        u = _acb_vec_init(n);
        v = _acb_vec_init(n);

        if (regularized == 2)
        {
            w = _acb_vec_init(n);
            acb_neg(t, s);
            _acb_poly_pow_acb_series(w, h, hlen, t, n, prec);
        }

        /* gamma(s, h(x)) = integral(h'(x) h(x)^(s-1) exp(-h(x)) */
        acb_sub_ui(u, s, 1, prec);
        _acb_poly_pow_acb_series(t, h, hlen, u, n, prec);
        _acb_poly_derivative(u, h, hlen, prec);
        _acb_poly_mullow(v, t, n, u, hlen - 1, n, prec);
        _acb_vec_neg(t, h, hlen);
        _acb_poly_exp_series(t, t, hlen, n, prec);
        _acb_poly_mullow(g, v, n, t, n, n, prec);
        _acb_poly_integral(g, g, n, prec);

        if (regularized == 1)
        {
            acb_rgamma(t, s, prec);
            _acb_vec_scalar_mul(g, g, n, t, prec);
        }
        else if (regularized == 2)
        {
            acb_rgamma(t, s, prec);
            _acb_vec_scalar_mul(g, g, n, t, prec);
            _acb_vec_set(u, g, n);
            _acb_poly_mullow(g, u, n, w, n, n, prec);
            _acb_vec_clear(w, n);
        }

        _acb_vec_clear(t, n);
        _acb_vec_clear(u, n);
        _acb_vec_clear(v, n);
    }

    acb_swap(g, c);
    acb_clear(c);
}


void
acb_hypgeom_gamma_lower_series(acb_poly_t g, const acb_t s,
        const acb_poly_t h, int regularized, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        acb_poly_zero(g);
        return;
    }

    acb_poly_fit_length(g, n);

    if (hlen == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_hypgeom_gamma_lower_series(g->coeffs, s, t, 1, regularized, n, prec);
        acb_clear(t);
    }
    else
    {
        _acb_hypgeom_gamma_lower_series(
                g->coeffs, s, h->coeffs, hlen, regularized, n, prec);
    }

    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}
