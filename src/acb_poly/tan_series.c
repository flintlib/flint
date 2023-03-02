/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

#define TAN_NEWTON_CUTOFF 20

void
_acb_poly_tan_series(acb_ptr g,
    acb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_tan(g, h, prec);
        _acb_vec_zero(g + 1, len - 1);
    }
    else if (len == 2)
    {
        acb_t t;
        acb_init(t);
        acb_tan(g, h, prec);
        acb_mul(t, g, g, prec);
        acb_add_ui(t, t, 1, prec);
        acb_mul(g + 1, t, h + 1, prec);  /* safe since hlen >= 2 */
        acb_clear(t);
    }
    else
    {
        acb_ptr t, u;

        t = _acb_vec_init(2 * len);
        u = t + len;

        NEWTON_INIT(TAN_NEWTON_CUTOFF, len)

        NEWTON_BASECASE(n)
        _acb_poly_sin_cos_series_basecase(t, u, h, hlen, n, prec, 0);
        _acb_poly_div_series(g, t, n, u, n, n, prec);
        NEWTON_END_BASECASE

        NEWTON_LOOP(m, n)
        _acb_poly_mullow(u, g, m, g, m, n, prec);
        acb_add_ui(u, u, 1, prec);
        _acb_poly_atan_series(t, g, m, n, prec);
        _acb_poly_sub(t + m, h + m, FLINT_MAX(0, hlen - m), t + m, n - m, prec);
        _acb_poly_mullow(g + m, u, n, t + m, n - m, n - m, prec);
        NEWTON_END_LOOP

        NEWTON_END

        _acb_vec_clear(t, 2 * len);
    }
}

void
acb_poly_tan_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
{
    if (h->length == 0 || n == 0)
    {
        acb_poly_zero(g);
        return;
    }

    if (g == h)
    {
        acb_poly_t t;
        acb_poly_init(t);
        acb_poly_tan_series(t, h, n, prec);
        acb_poly_swap(g, t);
        acb_poly_clear(t);
        return;
    }

    acb_poly_fit_length(g, n);
    _acb_poly_tan_series(g->coeffs, h->coeffs, h->length, n, prec);
    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}

