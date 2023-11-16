/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

slong arf_get_si(const arf_t x, arf_rnd_t rnd);

void _arb_poly_lgamma_series_at_one(arb_ptr u, slong len, slong prec);

void arb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const arb_t x, int use_reflect, int digamma, slong prec);

void _arb_poly_gamma_stirling_eval(arb_ptr res, const arb_t z, slong n, slong num, slong prec);


void
_arb_poly_rgamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec)
{
    int reflect, isint;
    slong i, rflen, r, n, wp;
    arb_ptr t, u, v;
    arb_struct f[2];

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        arb_rgamma(res, h, prec);
        _arb_vec_zero(res + 1, len - 1);
        return;
    }

    isint = arb_is_int(h);

    if (len <= 2 && isint && arf_sgn(arb_midref(h)) <= 0)
    {
        int even = arf_is_int_2exp_si(arb_midref(h), 1);

        /* use res[0] as tmp to allow aliasing */
        arb_sub_ui(res, h, 1, prec);
        arb_neg(res, res);
        arb_gamma(res, res, prec);
        arb_mul(res + 1, h + 1, res, prec);
        if (!even)
            arb_neg(res + 1, res + 1);

        arb_zero(res);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    t = _arb_vec_init(len);
    u = _arb_vec_init(len);
    v = _arb_vec_init(len);
    arb_init(f);
    arb_init(f + 1);

    /* use zeta values at small integers */
    if (isint && (arf_cmpabs_ui(arb_midref(h), prec / 2) < 0))
    {
        r = arf_get_si(arb_midref(h), ARF_RND_DOWN);

        _arb_poly_lgamma_series_at_one(u, len, wp);

        _arb_vec_neg(u, u, len);
        _arb_poly_exp_series(t, u, len, len, wp);

        if (r == 1)
        {
            _arb_vec_swap(v, t, len);
        }
        else if (r <= 0)
        {
            arb_set(f, h);
            arb_one(f + 1);
            rflen = FLINT_MIN(len, 2 - r);
            _arb_poly_rising_ui_series(u, f, FLINT_MIN(2, len), 1 - r, rflen, wp);
            _arb_poly_mullow(v, t, len, u, rflen, len, wp);
        }
        else
        {
            arb_one(f);
            arb_one(f + 1);
            rflen = FLINT_MIN(len, r);
            _arb_poly_rising_ui_series(v, f, FLINT_MIN(2, len), r - 1, rflen, wp);

            /* TODO: use div_series? */
            _arb_poly_inv_series(u, v, rflen, len, wp);
            _arb_poly_mullow(v, t, len, u, len, len, wp);
        }
    }
    else
    {
        /* otherwise use Stirling series */
        arb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, h, 1, 0, wp);

        /* rgamma(h) = (gamma(1-h+r) sin(pi h)) / (rf(1-h, r) * pi), h = h0 + t*/
        if (reflect)
        {
            /* u = gamma(r+1-h) */
            arb_sub_ui(f, h, r + 1, wp);
            arb_neg(f, f);
            _arb_poly_gamma_stirling_eval(t, f, n, len, wp);
            _arb_poly_exp_series(u, t, len, len, wp);
            for (i = 1; i < len; i += 2)
                arb_neg(u + i, u + i);

            /* v = sin(pi x) */
            arb_set(f, h);
            arb_one(f + 1);
            _arb_poly_sin_pi_series(v, f, 2, len, wp);

            _arb_poly_mullow(t, u, len, v, len, len, wp);

            /* rf(1-h,r) * pi */
            if (r == 0)
            {
                arb_const_pi(u, wp);
                _arb_vec_scalar_div(v, t, len, u, wp);
            }
            else
            {
                arb_sub_ui(f, h, 1, wp);
                arb_neg(f, f);
                arb_set_si(f + 1, -1);
                rflen = FLINT_MIN(len, r + 1);
                _arb_poly_rising_ui_series(v, f, FLINT_MIN(2, len), r, rflen, wp);
                arb_const_pi(u, wp);
                _arb_vec_scalar_mul(v, v, rflen, u, wp);

                /* divide by rising factorial */
                /* TODO: might better to use div_series, when it has a good basecase */
                _arb_poly_inv_series(u, v, rflen, len, wp);
                _arb_poly_mullow(v, t, len, u, len, len, wp);
            }
        }
        else
        {
            /* rgamma(h) = rgamma(h+r) rf(h,r) */
            if (r == 0)
            {
                arb_add_ui(f, h, r, wp);
                _arb_poly_gamma_stirling_eval(t, f, n, len, wp);
                _arb_vec_neg(t, t, len);
                _arb_poly_exp_series(v, t, len, len, wp);
            }
            else
            {
                arb_set(f, h);
                arb_one(f + 1);
                rflen = FLINT_MIN(len, r + 1);
                _arb_poly_rising_ui_series(t, f, FLINT_MIN(2, len), r, rflen, wp);

                arb_add_ui(f, h, r, wp);
                _arb_poly_gamma_stirling_eval(v, f, n, len, wp);
                _arb_vec_neg(v, v, len);
                _arb_poly_exp_series(u, v, len, len, wp);

                _arb_poly_mullow(v, u, len, t, rflen, len, wp);
            }
        }
    }

    /* compose with nonconstant part */
    arb_zero(t);
    _arb_vec_set(t + 1, h + 1, hlen - 1);
    _arb_poly_compose_series(res, v, len, t, hlen, len, prec);

    arb_clear(f);
    arb_clear(f + 1);
    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);
    _arb_vec_clear(v, len);
}

void
arb_poly_rgamma_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
{
    if (f->length == 0 || n == 0)
    {
        arb_poly_zero(res);
    }
    else
    {
        arb_poly_fit_length(res, n);
        _arb_poly_rgamma_series(res->coeffs, f->coeffs, f->length, n, prec);
        _arb_poly_set_length(res, n);
        _arb_poly_normalise(res);
    }
}

