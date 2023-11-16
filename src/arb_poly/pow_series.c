/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_pow_series(arb_ptr h,
    arb_srcptr f, slong flen,
    arb_srcptr g, slong glen, slong len, slong prec)
{
    if (glen == 1)
    {
        _arb_poly_pow_arb_series(h, f, flen, g, len, prec);
        return;
    }

    /* f^g = exp(g * log(f)) */
    if (flen == 1)
    {
        arb_t t;
        arb_init(t);
        arb_log(t, f, prec);
        _arb_vec_scalar_mul(h, g, glen, t, prec);
        _arb_poly_exp_series(h, h, glen, len, prec);
        arb_clear(t);
    }
    else
    {
        arb_ptr t;
        t = _arb_vec_init(len);
        _arb_poly_log_series(t, f, flen, len, prec);
        _arb_poly_mullow(h, t, len, g, glen, len, prec);
        _arb_poly_exp_series(h, h, len, len, prec);
        _arb_vec_clear(t, len);

    }
}

void
arb_poly_pow_series(arb_poly_t h,
    const arb_poly_t f, const arb_poly_t g, slong len, slong prec)
{
    slong flen, glen;

    flen = f->length;
    glen = g->length;

    flen = FLINT_MIN(flen, len);
    glen = FLINT_MIN(glen, len);

    if (len == 0)
    {
        arb_poly_zero(h);
        return;
    }

    if (glen == 0)
    {
        arb_poly_one(h);
        return;
    }

    if (flen == 0)
    {
        arb_poly_zero(h);
        return;
    }

    if (flen == 1 && glen == 1)
    {
        arb_poly_fit_length(h, 1);
        arb_pow(h->coeffs, f->coeffs, g->coeffs, prec);
        _arb_poly_set_length(h, 1);
        _arb_poly_normalise(h);
        return;
    }

    if (f == h || g == h)
    {
        arb_poly_t t;
        arb_poly_init2(t, len);
        _arb_poly_pow_series(t->coeffs,
            f->coeffs, flen, g->coeffs, glen, len, prec);
        _arb_poly_set_length(t, len);
        _arb_poly_normalise(t);
        arb_poly_swap(t, h);
        arb_poly_clear(t);
    }
    else
    {
        arb_poly_fit_length(h, len);
        _arb_poly_pow_series(h->coeffs,
            f->coeffs, flen, g->coeffs, glen, len, prec);
        _arb_poly_set_length(h, len);
        _arb_poly_normalise(h);
    }
}

