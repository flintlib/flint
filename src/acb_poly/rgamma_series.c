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

void acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t x, int use_reflect, int digamma, slong prec);

void
_acb_poly_gamma_stirling_eval(acb_ptr res, const acb_t z, slong n, slong num, slong prec);

void
_acb_poly_rgamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong len, slong prec)
{
    int reflect;
    slong i, rflen, r, n, wp;
    acb_ptr t, u, v;
    acb_struct f[2];

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_rgamma(res, h, prec);
        _acb_vec_zero(res + 1, len - 1);
        return;
    }

    /* use real code for real input */
    if (_acb_vec_is_real(h, hlen))
    {
        arb_ptr tmp = _arb_vec_init(len);
        for (i = 0; i < hlen; i++)
            arb_set(tmp + i, acb_realref(h + i));
        _arb_poly_rgamma_series(tmp, tmp, hlen, len, prec);
        for (i = 0; i < len; i++)
            acb_set_arb(res + i, tmp + i);
        _arb_vec_clear(tmp, len);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);
    v = _acb_vec_init(len);
    acb_init(f);
    acb_init(f + 1);

    /* otherwise use Stirling series */
    acb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, h, 1, 0, wp);

    /* rgamma(h) = (gamma(1-h+r) sin(pi h)) / (rf(1-h, r) * pi), h = h0 + t*/
    if (reflect)
    {
        /* u = gamma(r+1-h) */
        acb_sub_ui(f, h, r + 1, wp);
        acb_neg(f, f);
        _acb_poly_gamma_stirling_eval(t, f, n, len, wp);
        _acb_poly_exp_series(u, t, len, len, wp);
        for (i = 1; i < len; i += 2)
            acb_neg(u + i, u + i);

        /* v = sin(pi x) */
        acb_set(f, h);
        acb_one(f + 1);
        _acb_poly_sin_pi_series(v, f, 2, len, wp);

        _acb_poly_mullow(t, u, len, v, len, len, wp);

        /* rf(1-h,r) * pi */
        if (r == 0)
        {
            acb_const_pi(u, wp);
            _acb_vec_scalar_div(v, t, len, u, wp);
        }
        else
        {
            acb_sub_ui(f, h, 1, wp);
            acb_neg(f, f);
            acb_set_si(f + 1, -1);
            rflen = FLINT_MIN(len, r + 1);
            _acb_poly_rising_ui_series(v, f, FLINT_MIN(2, len), r, rflen, wp);
            acb_const_pi(u, wp);
            _acb_vec_scalar_mul(v, v, rflen, u, wp);

            /* divide by rising factorial */
            /* TODO: might better to use div_series, when it has a good basecase */
            _acb_poly_inv_series(u, v, rflen, len, wp);
            _acb_poly_mullow(v, t, len, u, len, len, wp);
        }
    }
    else
    {
        /* rgamma(h) = rgamma(h+r) rf(h,r) */
        if (r == 0)
        {
            acb_add_ui(f, h, r, wp);
            _acb_poly_gamma_stirling_eval(t, f, n, len, wp);
            _acb_vec_neg(t, t, len);
            _acb_poly_exp_series(v, t, len, len, wp);
        }
        else
        {
            acb_set(f, h);
            acb_one(f + 1);
            rflen = FLINT_MIN(len, r + 1);
            _acb_poly_rising_ui_series(t, f, FLINT_MIN(2, len), r, rflen, wp);

            acb_add_ui(f, h, r, wp);
            _acb_poly_gamma_stirling_eval(v, f, n, len, wp);
            _acb_vec_neg(v, v, len);
            _acb_poly_exp_series(u, v, len, len, wp);

            _acb_poly_mullow(v, u, len, t, rflen, len, wp);
        }
    }

    /* compose with nonconstant part */
    acb_zero(t);
    _acb_vec_set(t + 1, h + 1, hlen - 1);
    _acb_poly_compose_series(res, v, len, t, hlen, len, prec);

    acb_clear(f);
    acb_clear(f + 1);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
    _acb_vec_clear(v, len);
}

void
acb_poly_rgamma_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
{
    if (f->length == 0 || n == 0)
    {
        acb_poly_zero(res);
    }
    else
    {
        acb_poly_fit_length(res, n);
        _acb_poly_rgamma_series(res->coeffs, f->coeffs, f->length, n, prec);
        _acb_poly_set_length(res, n);
        _acb_poly_normalise(res);
    }
}

