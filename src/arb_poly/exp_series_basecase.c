/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

static void
_arb_poly_exp_series_basecase_rec(arb_ptr f, arb_ptr a,
        arb_srcptr h, slong hlen, slong n, slong prec)
{
    slong k;

    arb_t s;
    arb_init(s);

    arb_exp(f, h, prec);

    for (k = 1; k < hlen; k++)
        arb_mul_ui(a + k, h + k, k, prec);

    for (k = 1; k < n; k++)
    {
        arb_dot(s, NULL, 0, a + 1, 1, f + k - 1, -1, FLINT_MIN(k, hlen - 1), prec);
        arb_div_ui(f + k, s, k, prec);
    }

    arb_clear(s);
}

void
_arb_poly_exp_series_basecase(arb_ptr f,
        arb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(n, hlen);

    if (n < 20 || hlen < 0.9 * n || prec <= 2 * FLINT_BITS || n < 1000.0 / log(prec + 10) - 70)
    {
        arb_ptr t = _arb_vec_init(hlen);
        _arb_poly_exp_series_basecase_rec(f, t, h, hlen, n, prec);
        _arb_vec_clear(t, hlen);
    }
    else
    {
        slong m, v;
        arb_ptr t, u;

        m = (n + 2) / 3;
        v = m * 2;

        t = _arb_vec_init(n);
        u = _arb_vec_init(n - m);

        _arb_poly_mullow(t, h + m, hlen - m, h + m, hlen - m, n - v, prec);
        _arb_vec_scalar_mul_2exp_si(t, t, n - v, -1);

        _arb_vec_set(u, h + m, v - m);
        _arb_poly_add(u + v - m, t, n - v, h + v, hlen - v, prec);
        _arb_poly_exp_series_basecase_rec(f, t, h, m, n, prec);
        _arb_poly_mullow(t, f, n, u, n - m, n - m, prec);
        _arb_poly_add(f + m, f + m, n - m, t, n - m, prec);

        _arb_vec_clear(t, n);
        _arb_vec_clear(u, n - m);
    }
}

void
arb_poly_exp_series_basecase(arb_poly_t f, const arb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        arb_poly_zero(f);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_one(f);
        return;
    }

    arb_poly_fit_length(f, n);
    _arb_poly_exp_series_basecase(f->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(f, n);
    _arb_poly_normalise(f);
}
