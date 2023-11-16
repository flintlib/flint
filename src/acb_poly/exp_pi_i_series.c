/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_exp_pi_i_series(acb_ptr f, acb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_exp_pi_i(f, h, prec);
        _acb_vec_zero(f + 1, len - 1);
    }
    else if (len == 2)
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, prec);
        acb_exp_pi_i(f, h, prec);
        acb_mul_arb(f + 1, h + 1, t, prec);
        acb_mul_onei(f + 1, f + 1);
        acb_mul(f + 1, f + 1, f + 0, prec);
        arb_clear(t);
    }
    else
    {
        acb_ptr t;
        t = _acb_vec_init(hlen + 1);
        acb_const_pi(t, prec);
        acb_mul_onei(t, t);
        _acb_vec_scalar_mul(t + 1, h + 1, hlen - 1, t, prec);
        acb_zero(t);
        acb_exp_pi_i(t + hlen, h, prec);
        _acb_poly_exp_series(f, t, hlen, len, prec);
        _acb_vec_scalar_mul(f, f, len, t + hlen, prec);
        _acb_vec_clear(t, hlen + 1);
    }
}

void
acb_poly_exp_pi_i_series(acb_poly_t f, const acb_poly_t h, slong n, slong prec)
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

    if (hlen == 1)
        n = 1;

    acb_poly_fit_length(f, n);
    _acb_poly_exp_pi_i_series(f->coeffs, h->coeffs, hlen, n, prec);
    _acb_poly_set_length(f, n);
    _acb_poly_normalise(f);
}

