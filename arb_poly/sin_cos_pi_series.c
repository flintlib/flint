/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        arb_sin_cos_pi(s, c, h, prec);
        _arb_vec_zero(s + 1, n - 1);
        _arb_vec_zero(c + 1, n - 1);
    }
    else if (n == 2)
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, prec);
        arb_mul(t, t, h + 1, prec);
        arb_sin_cos_pi(s, c, h, prec);
        arb_mul(s + 1, c, t, prec);
        arb_neg(t, t);
        arb_mul(c + 1, s, t, prec);
        arb_clear(t);
    }
    else
    {
        slong cutoff;

        if (prec <= 128)
        {
            cutoff = 1400;
        }
        else
        {
            cutoff = 100000 / pow(log(prec), 3);
            cutoff = FLINT_MIN(cutoff, 700);
        }

        if (hlen < cutoff)
            _arb_poly_sin_cos_series_basecase(s, c, h, hlen, n, prec, 1);
        else
            _arb_poly_sin_cos_series_tangent(s, c, h, hlen, n, prec, 1);
    }
}

void
arb_poly_sin_cos_pi_series(arb_poly_t s, arb_poly_t c,
                                    const arb_poly_t h, slong n, slong prec)
{
    slong hlen = h->length;

    if (n == 0)
    {
        arb_poly_zero(s);
        arb_poly_zero(c);
        return;
    }

    if (hlen == 0)
    {
        arb_poly_zero(s);
        arb_poly_one(c);
        return;
    }

    if (hlen == 1)
        n = 1;

    arb_poly_fit_length(s, n);
    arb_poly_fit_length(c, n);
    _arb_poly_sin_cos_pi_series(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(s, n);
    _arb_poly_normalise(s);
    _arb_poly_set_length(c, n);
    _arb_poly_normalise(c);
}

