/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_sqrt_series(acb_ptr g,
    acb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    while (hlen > 0 && acb_is_zero(h + hlen - 1))
        hlen--;

    if (hlen <= 1)
    {
        acb_sqrt(g, h, prec);
        _acb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        acb_sqrt(g, h, prec);
        acb_div(g + 1, h + 1, h, prec);
        acb_mul(g + 1, g + 1, g, prec);
        acb_mul_2exp_si(g + 1, g + 1, -1);
    }
    else if (_acb_vec_is_zero(h + 1, hlen - 2))
    {
        acb_t t;
        acb_init(t);
        arf_set_si_2exp_si(arb_midref(acb_realref(t)), 1, -1);
        _acb_poly_binomial_pow_acb_series(g, h, hlen, t, len, prec);
        acb_clear(t);
    }
    else
    {
        /* Karp-Markstein: with g = h^(-1/2) + O(x^m) and v = h g mod x^m,
           sqrt(h) = v + g (h - v^2) / 2 + O(x^(2m)). */
        acb_ptr t, u;
        slong m, tlen, r1, r2, rlen;

        m = (len + 1) / 2;
        t = _acb_vec_init(m + len);
        u = t + m;

        _acb_poly_rsqrt_series(t, h, hlen, m, prec);
        _acb_poly_mullow(g, t, m, h, FLINT_MIN(hlen, m), m, prec);

        tlen = FLINT_MIN(2 * m - 1, len);
        r1 = FLINT_MAX(0, FLINT_MIN(hlen - m, len - m));
        r2 = FLINT_MAX(0, tlen - m);
        rlen = FLINT_MAX(r1, r2);

        if (r2 > 0)
            _acb_poly_mulmid(u + m, g, m, g, m, m, tlen, prec);
        _acb_poly_sub(u, h + m, r1, u + m, r2, prec);

        _acb_poly_mullow(g + m, t, len - m, u, rlen, len - m, prec);
        _acb_vec_scalar_mul_2exp_si(g + m, g + m, len - m, -1);

        _acb_vec_clear(t, m + len);
    }
}

void
acb_poly_sqrt_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(g);
        return;
    }

    if (g == h)
    {
        acb_poly_t t;
        acb_poly_init(t);
        acb_poly_sqrt_series(t, h, n, prec);
        acb_poly_swap(g, t);
        acb_poly_clear(t);
        return;
    }

    acb_poly_fit_length(g, n);
    if (h->length == 0)
        _acb_vec_indeterminate(g->coeffs, n);
    else
        _acb_poly_sqrt_series(g->coeffs, h->coeffs, h->length, n, prec);
    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}
