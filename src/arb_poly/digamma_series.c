/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void arb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const arb_t x, int use_reflect, int digamma, slong prec);

void _arb_poly_gamma_stirling_eval2(arb_ptr res, const arb_t z, slong n, slong num, int diff, slong prec);

void
_arb_poly_digamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec)
{
    int reflect;
    slong i, r, n, rflen, wp;
    arb_t zr;
    arb_ptr t, u, v;

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        arb_digamma(res, h, prec);
        if (arb_is_finite(res))
            _arb_vec_zero(res + 1, len - 1);
        else
            _arb_vec_indeterminate(res + 1, len - 1);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(prec);

    t = _arb_vec_init(len + 1);
    u = _arb_vec_init(len + 1);
    v = _arb_vec_init(len + 1);
    arb_init(zr);

    /* use zeta values at small integers */
    if (arb_is_int(h) && (arf_cmpabs_ui(arb_midref(h), prec / 2) < 0))
    {
        r = arf_get_si(arb_midref(h), ARF_RND_DOWN);

        if (r <= 0)
        {
            _arb_vec_indeterminate(res, len);
            goto cleanup;
        }
        else if (r == 1)
        {
            arb_const_euler(u, prec);
            arb_zeta_ui_vec(u + 1, 2, len - 1, prec);
            for (i = 0; i < len; i += 2)
                arb_neg(u + i, u + i);
        }
        else
        {
            arb_one(v);
            arb_one(v + 1);
            rflen = FLINT_MIN(len + 1, r);
            _arb_poly_rising_ui_series(u, v, 2, r - 1, rflen, wp);
            _arb_poly_derivative(v, u, rflen, wp);
            _arb_poly_div_series(t, v, rflen - 1, u, rflen, len, wp);

            arb_const_euler(u, prec);
            arb_zeta_ui_vec(u + 1, 2, len - 1, prec);
            for (i = 0; i < len; i += 2)
                arb_neg(u + i, u + i);

            _arb_vec_add(u, u, t, len, wp);
        }
    }
    else
    {
        /* use Stirling series */
        arb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, h, 1, 1, wp);

        /* psi(x) = psi((1-x)+r) - h(1-x,r) - pi*cot(pi*x) */
        if (reflect)
        {
            if (r != 0) /* otherwise t = 0 */
            {
                arb_sub_ui(v, h, 1, wp);
                arb_neg(v, v);
                arb_one(v + 1);
                rflen = FLINT_MIN(len + 1, r + 1);
                _arb_poly_rising_ui_series(u, v, 2, r, rflen, wp);
                _arb_poly_derivative(v, u, rflen, wp);
                _arb_poly_div_series(t, v, rflen - 1, u, rflen, len, wp);
                for (i = 1; i < len; i += 2)
                    arb_neg(t + i, t + i);
            }

            arb_sub_ui(zr, h, r + 1, wp);
            arb_neg(zr, zr);
            _arb_poly_gamma_stirling_eval2(u, zr, n, len + 1, 1, wp);
            for (i = 1; i < len; i += 2)
                arb_neg(u + i, u + i);

            _arb_vec_sub(u, u, t, len, wp);

            arb_set(t, h);
            arb_one(t + 1);
            _arb_poly_cot_pi_series(t, t, 2, len, wp);
            arb_const_pi(v, wp);
            _arb_vec_scalar_mul(t, t, len, v, wp);

            _arb_vec_sub(u, u, t, len, wp);
        }
        else
        {
            if (r == 0)
            {
                arb_add_ui(zr, h, r, wp);
                _arb_poly_gamma_stirling_eval2(u, zr, n, len + 1, 1, wp);
            }
            else
            {
                arb_set(v, h);
                arb_one(v + 1);
                rflen = FLINT_MIN(len + 1, r + 1);
                _arb_poly_rising_ui_series(u, v, 2, r, rflen, wp);
                _arb_poly_derivative(v, u, rflen, wp);
                _arb_poly_div_series(t, v, rflen - 1, u, rflen, len, wp);

                arb_add_ui(zr, h, r, wp);
                _arb_poly_gamma_stirling_eval2(u, zr, n, len + 1, 1, wp);

                _arb_vec_sub(u, u, t, len, wp);
            }
        }
    }

    /* compose with nonconstant part */
    arb_zero(t);
    _arb_vec_set(t + 1, h + 1, hlen - 1);
    _arb_poly_compose_series(res, u, len, t, hlen, len, prec);

cleanup:
    arb_clear(zr);
    _arb_vec_clear(t, len + 1);
    _arb_vec_clear(u, len + 1);
    _arb_vec_clear(v, len + 1);
}

void
arb_poly_digamma_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
{
    arb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _arb_vec_indeterminate(res->coeffs, n);
    else
        _arb_poly_digamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

