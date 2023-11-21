/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

static void
_acb_poly_exp_series_basecase_rec(acb_ptr f, acb_ptr a,
        acb_srcptr h, slong hlen, slong n, slong prec)
{
    slong k;

    acb_t s;
    acb_init(s);

    acb_exp(f, h, prec);

    for (k = 1; k < hlen; k++)
        acb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        acb_dot(s, NULL, 0, a + 1, 1, f + k - 1, -1, FLINT_MIN(k, hlen - 1), prec);
        acb_div_ui(f + k, s, k, prec);
    }

    acb_clear(s);
}

void
_acb_poly_exp_series_basecase(acb_ptr f,
        acb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(n, hlen);

    if (n < 20 || hlen < 0.9 * n || prec <= 2 * FLINT_BITS || n < 1000.0 / log(prec + 10) - 70)
    {
        acb_ptr t = _acb_vec_init(hlen);
        _acb_poly_exp_series_basecase_rec(f, t, h, hlen, n, prec);
        _acb_vec_clear(t, hlen);
    }
    else
    {
        slong m, v;
        acb_ptr t, u;

        m = (n + 2) / 3;
        v = m * 2;

        t = _acb_vec_init(n);
        u = _acb_vec_init(n - m);

        _acb_poly_mullow(t, h + m, hlen - m, h + m, hlen - m, n - v, prec);
        _acb_vec_scalar_mul_2exp_si(t, t, n - v, -1);

        _acb_vec_set(u, h + m, v - m);
        _acb_poly_add(u + v - m, t, n - v, h + v, hlen - v, prec);
        _acb_poly_exp_series_basecase_rec(f, t, h, m, n, prec);
        _acb_poly_mullow(t, f, n, u, n - m, n - m, prec);
        _acb_poly_add(f + m, f + m, n - m, t, n - m, prec);

        _acb_vec_clear(t, n);
        _acb_vec_clear(u, n - m);
    }
}

void
acb_poly_exp_series_basecase(acb_poly_t f, const acb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        acb_poly_zero(f);
        return;
    }

    if (hlen == 0)
    {
        acb_poly_one(f);
        return;
    }

    acb_poly_fit_length(f, n);
    _acb_poly_exp_series_basecase(f->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(f, n);
    _acb_poly_normalise(f);
}
