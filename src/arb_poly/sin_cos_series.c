/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

#ifdef __GNUC__
# define log __builtin_log
# define pow __builtin_pow
#else
# include <math.h>
#endif

static void
__arb_poly_sin_cos_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, int times_pi, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        if (times_pi)
            arb_sin_cos_pi(s, c, h, prec);
        else
            arb_sin_cos(s, c, h, prec);

        _arb_vec_zero(s + 1, n - 1);
        _arb_vec_zero(c + 1, n - 1);
    }
    else if (n == 2)
    {
        arb_t t;
        arb_init(t);
        if (times_pi)
        {
            arb_const_pi(t, prec);
            arb_mul(t, t, h + 1, prec);
            arb_sin_cos_pi(s, c, h, prec);
        }
        else
        {
            arb_set(t, h + 1);
            arb_sin_cos(s, c, h, prec);
        }

        arb_mul(s + 1, c, t, prec);
        arb_neg(t, t);
        arb_mul(c + 1, s, t, prec);
        arb_clear(t);
    }
    else
    {
        slong cutoff;
        gr_ctx_t ctx;
        int status;

        if (prec <= 128)
        {
            cutoff = 1400;
        }
        else
        {
            cutoff = 100000 / pow(log(prec), 3);
            cutoff = FLINT_MIN(cutoff, 700);
        }

        gr_ctx_init_real_arb(ctx, prec);

        if (hlen < cutoff)
            status = _gr_poly_sin_cos_series_basecase(s, c, h, hlen, n, times_pi, ctx);
        else
            status = _gr_poly_sin_cos_series_tangent(s, c, h, hlen, n, times_pi, ctx);

        if (status != GR_SUCCESS)
        {
            _arb_vec_indeterminate(s, n);
            _arb_vec_indeterminate(c, n);
        }
    }
}

void
_arb_poly_sin_cos_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
{
    __arb_poly_sin_cos_series(s, c, h, hlen, n, 0, prec);
}

void
_arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
{
    __arb_poly_sin_cos_series(s, c, h, hlen, n, 1, prec);
}

void
arb_poly_sin_cos_series(arb_poly_t s, arb_poly_t c,
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
    _arb_poly_sin_cos_series(s->coeffs, c->coeffs, h->coeffs, hlen, n, prec);
    _arb_poly_set_length(s, n);
    _arb_poly_normalise(s);
    _arb_poly_set_length(c, n);
    _arb_poly_normalise(c);
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
