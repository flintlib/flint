/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

#ifdef __GNUC__
# define log __builtin_log
# define pow __builtin_pow
#else
# include <math.h>
#endif

void
_acb_poly_sinh_cosh_series(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        acb_sinh_cosh(s, c, h, prec);
        _acb_vec_zero(s + 1, n - 1);
        _acb_vec_zero(c + 1, n - 1);
    }
    else if (n == 2)
    {
        acb_t t;
        acb_init(t);
        acb_set(t, h + 1);
        acb_sinh_cosh(s, c, h, prec);
        acb_mul(s + 1, c, t, prec);
        acb_mul(c + 1, s, t, prec);
        acb_clear(t);
    }
    else
    {
        slong cutoff;

        if (prec <= 128)
            cutoff = 400;
        else
            cutoff = 30000 / pow(log(prec), 3);

        if (hlen < cutoff)
            _acb_poly_sinh_cosh_series_basecase(s, c, h, hlen, n, prec);
        else
            _acb_poly_sinh_cosh_series_exponential(s, c, h, hlen, n, prec);
    }
}

void
acb_poly_sinh_cosh_series(acb_poly_t s, acb_poly_t c,
                                    const acb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        acb_poly_zero(s);
        acb_poly_zero(c);
        return;
    }

    if (hlen == 0)
    {
        acb_poly_zero(s);
        acb_poly_one(c);
        return;
    }

    if (hlen == 1)
        n = 1;

    acb_poly_fit_length(s, n);
    acb_poly_fit_length(c, n);
    _acb_poly_sinh_cosh_series(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(s, n);
    _acb_poly_normalise(s);
    _acb_poly_set_length(c, n);
    _acb_poly_normalise(c);
}
