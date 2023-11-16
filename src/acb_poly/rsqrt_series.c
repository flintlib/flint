/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"

void
_acb_poly_rsqrt_series(acb_ptr g,
    acb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    while (hlen > 0 && acb_is_zero(h + hlen - 1))
        hlen--;

    if (hlen <= 1)
    {
        acb_rsqrt(g, h, prec);
        _acb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        acb_rsqrt(g, h, prec);
        acb_div(g + 1, h + 1, h, prec);
        acb_mul(g + 1, g + 1, g, prec);
        acb_mul_2exp_si(g + 1, g + 1, -1);
        acb_neg(g + 1, g + 1);
    }
    else if (_acb_vec_is_zero(h + 1, hlen - 2))
    {
        acb_t t;
        acb_init(t);
        arf_set_si_2exp_si(arb_midref(acb_realref(t)), -1, -1);
        _acb_poly_binomial_pow_acb_series(g, h, hlen, t, len, prec);
        acb_clear(t);
    }
    else
    {
        acb_ptr t, u;
        slong tlen;
        t = _acb_vec_init(2 * len);
        u = t + len;

        acb_rsqrt(g, h, prec);

        NEWTON_INIT(1, len)

        NEWTON_LOOP(m, n)
        tlen = FLINT_MIN(2 * m - 1, n);
        _acb_poly_mullow(t, g, m, g, m, tlen, prec);
        _acb_poly_mullow(u, g, m, t, tlen, n, prec);
        _acb_poly_mullow(t, u, n, h, hlen, n, prec);
        _acb_vec_scalar_mul_2exp_si(g + m, t + m, n - m, -1);
        _acb_vec_neg(g + m, g + m, n - m);
        NEWTON_END_LOOP

        NEWTON_END

        _acb_vec_clear(t, 2 * len);
    }
}

void
acb_poly_rsqrt_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
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
        acb_poly_rsqrt_series(t, h, n, prec);
        acb_poly_swap(g, t);
        acb_poly_clear(t);
        return;
    }

    acb_poly_fit_length(g, n);
    if (h->length == 0)
        _acb_vec_indeterminate(g->coeffs, n);
    else
        _acb_poly_rsqrt_series(g->coeffs, h->coeffs, h->length, n, prec);
    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}

