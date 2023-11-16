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
#include "acb_hypgeom.h"

void
_acb_log_rising_correct_branch(acb_t t,
        const acb_t t_wrong, const acb_t z, ulong r, slong prec);

void acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t x, int use_reflect, int digamma, slong prec);

void
_acb_poly_gamma_stirling_eval(acb_ptr res, const acb_t z, slong n, slong num, slong prec);

void
_acb_poly_lgamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong len, slong prec)
{
    int reflect;
    slong i, r, n, wp;
    acb_t zr;
    acb_ptr t, u;

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_lgamma(res, h, prec);
        if (acb_is_finite(res))
            _acb_vec_zero(res + 1, len - 1);
        else
            _acb_vec_indeterminate(res + 1, len - 1);
        return;
    }

    if (len == 2)
    {
        acb_t v;
        acb_init(v);
        acb_set(v, h + 1);
        acb_digamma(res + 1, h, prec);
        acb_lgamma(res, h, prec);
        acb_mul(res + 1, res + 1, v, prec);
        acb_clear(v);
        return;
    }

    /* use real code for real input and output */
    if (_acb_vec_is_real(h, hlen) && arb_is_positive(acb_realref(h)))
    {
        arb_ptr tmp = _arb_vec_init(len);
        for (i = 0; i < hlen; i++)
            arb_set(tmp + i, acb_realref(h + i));
        _arb_poly_lgamma_series(tmp, tmp, hlen, len, prec);
        for (i = 0; i < len; i++)
            acb_set_arb(res + i, tmp + i);
        _arb_vec_clear(tmp, len);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);
    acb_init(zr);

    /* use Stirling series */
    acb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, h, 1, 0, wp);

    if (reflect)
    {
        /* log gamma(h+x) = log rf(1-(h+x), r) - log gamma(1-(h+x)+r) - log sin(pi (h+x)) + log(pi) */
        if (r != 0) /* otherwise t = 0 */
        {
            acb_sub_ui(u, h, 1, wp);
            acb_neg(u, u);
            acb_hypgeom_log_rising_ui_jet(t, u, r, len, wp);
            for (i = 1; i < len; i += 2)
                acb_neg(t + i, t + i);
        }

        acb_sub_ui(u, h, 1, wp);
        acb_neg(u, u);
        acb_add_ui(zr, u, r, wp);
        _acb_poly_gamma_stirling_eval(u, zr, n, len, wp);
        for (i = 1; i < len; i += 2)
            acb_neg(u + i, u + i);

        _acb_vec_sub(t, t, u, len, wp);

        /* log(sin) is unstable with large imaginary parts;
           cot_pi is implemented in a numerically stable way */
        acb_set(u, h);
        acb_one(u + 1);
        _acb_poly_cot_pi_series(u, u, 2, len - 1, wp);
        _acb_poly_integral(u, u, len, wp);
        acb_const_pi(u, wp);
        _acb_vec_scalar_mul(u + 1, u + 1, len - 1, u, wp);
        acb_log_sin_pi(u, h, wp);

        _acb_vec_sub(u, t, u, len, wp);

        acb_const_pi(t, wp); /* todo: constant for log pi */
        acb_log(t, t, wp);
        acb_add(u, u, t, wp);
    }
    else
    {
        /* log gamma(x) = log gamma(x+r) - log rf(x,r) */

        acb_add_ui(zr, h, r, wp);
        _acb_poly_gamma_stirling_eval(u, zr, n, len, wp);

        if (r != 0)
        {
            acb_hypgeom_log_rising_ui_jet(t, h, r, len, wp);
            _acb_vec_sub(u, u, t, len, wp);
        }
    }

    /* compose with nonconstant part */
    acb_zero(t);
    _acb_vec_set(t + 1, h + 1, hlen - 1);
    _acb_poly_compose_series(res, u, len, t, hlen, len, prec);

    acb_clear(zr);
    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
}

void
acb_poly_lgamma_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
{
    acb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _acb_vec_indeterminate(res->coeffs, n);
    else
        _acb_poly_lgamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

