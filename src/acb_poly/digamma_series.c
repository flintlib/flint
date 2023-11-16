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

void _acb_poly_gamma_stirling_eval2(acb_ptr res, const acb_t z, slong n, slong num, int diff, slong prec);

void
_acb_poly_digamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong len, slong prec)
{
    int reflect;
    slong i, r, n, rflen, wp;
    acb_t zr;
    acb_ptr t, u, v;

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_digamma(res, h, prec);
        if (acb_is_finite(res))
            _acb_vec_zero(res + 1, len - 1);
        else
            _acb_vec_indeterminate(res + 1, len - 1);
        return;
    }

    /* use real code for real input */
    if (_acb_vec_is_real(h, hlen))
    {
        arb_ptr tmp = _arb_vec_init(len);
        for (i = 0; i < hlen; i++)
            arb_set(tmp + i, acb_realref(h + i));
        _arb_poly_digamma_series(tmp, tmp, hlen, len, prec);
        for (i = 0; i < len; i++)
            acb_set_arb(res + i, tmp + i);
        _arb_vec_clear(tmp, len);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    t = _acb_vec_init(len + 1);
    u = _acb_vec_init(len + 1);
    v = _acb_vec_init(len + 1);
    acb_init(zr);

    /* use Stirling series */
    acb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, h, 1, 1, wp);

    /* psi(x) = psi((1-x)+r) - h(1-x,r) - pi*cot(pi*x) */
    if (reflect)
    {
        if (r != 0) /* otherwise t = 0 */
        {
            acb_sub_ui(v, h, 1, wp);
            acb_neg(v, v);
            acb_one(v + 1);
            rflen = FLINT_MIN(len + 1, r + 1);
            _acb_poly_rising_ui_series(u, v, 2, r, rflen, wp);
            _acb_poly_derivative(v, u, rflen, wp);
            _acb_poly_div_series(t, v, rflen - 1, u, rflen, len, wp);
            for (i = 1; i < len; i += 2)
                acb_neg(t + i, t + i);
        }

        acb_sub_ui(zr, h, r + 1, wp);
        acb_neg(zr, zr);
        _acb_poly_gamma_stirling_eval2(u, zr, n, len + 1, 1, wp);
        for (i = 1; i < len; i += 2)
            acb_neg(u + i, u + i);

        _acb_vec_sub(u, u, t, len, wp);

        acb_set(t, h);
        acb_one(t + 1);
        _acb_poly_cot_pi_series(t, t, 2, len, wp);
        acb_const_pi(v, wp);
        _acb_vec_scalar_mul(t, t, len, v, wp);

        _acb_vec_sub(u, u, t, len, wp);
    }
    else
    {
        if (r == 0)
        {
            acb_add_ui(zr, h, r, wp);
            _acb_poly_gamma_stirling_eval2(u, zr, n, len + 1, 1, wp);
        }
        else
        {
            acb_set(v, h);
            acb_one(v + 1);
            rflen = FLINT_MIN(len + 1, r + 1);
            _acb_poly_rising_ui_series(u, v, 2, r, rflen, wp);
            _acb_poly_derivative(v, u, rflen, wp);
            _acb_poly_div_series(t, v, rflen - 1, u, rflen, len, wp);

            acb_add_ui(zr, h, r, wp);
            _acb_poly_gamma_stirling_eval2(u, zr, n, len + 1, 1, wp);

            _acb_vec_sub(u, u, t, len, wp);
        }
    }

    /* compose with nonconstant part */
    acb_zero(t);
    _acb_vec_set(t + 1, h + 1, hlen - 1);
    _acb_poly_compose_series(res, u, len, t, hlen, len, prec);

    acb_clear(zr);
    _acb_vec_clear(t, len + 1);
    _acb_vec_clear(u, len + 1);
    _acb_vec_clear(v, len + 1);
}

void
acb_poly_digamma_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
{
    acb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _acb_vec_indeterminate(res->coeffs, n);
    else
        _acb_poly_digamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

