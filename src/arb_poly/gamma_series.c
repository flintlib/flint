/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void arb_gamma_stirling_bound(mag_ptr err, const arb_t x, slong k0, slong knum, slong n);

void arb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const arb_t x, int use_reflect, int digamma, slong prec);

void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);

void
_arb_poly_lgamma_series_at_one(arb_ptr u, slong len, slong prec)
{
    slong i;

    if (len > 0) arb_zero(u);
    if (len > 1) arb_const_euler(u + 1, prec);
    if (len > 2) arb_zeta_ui_vec(u + 2, 2, len - 2, prec);


    for (i = 2; i < len; i++)
        arb_div_ui(u + i, u + i, i, prec);

    for (i = 1; i < len; i += 2)
        arb_neg(u + i, u + i);
}

static void
bsplit(arb_ptr Q, arb_ptr T, const arb_t z, slong a, slong b, slong num, slong prec)
{
    if (b - a == 1)
    {
        arb_gamma_stirling_coeff(T, a, 0, prec);

        if (a == 1)
        {   /* (z + t) */
            arb_set(Q, z);
            if (num > 1) arb_one(Q + 1);
            if (num > 2) arb_zero(Q + 2);
        }
        else
        {   /* (z + t)^2 */
            arb_mul(Q, z, z, prec);  /* TODO: precompute */
            if (num > 1) arb_mul_2exp_si(Q + 1, z, 1);
            if (num > 2) arb_one(Q + 2);
        }
    }
    else
    {
        slong m, n1, n2, q1len, q2len, t1len, t2len, qlen, tlen, alloc;
        arb_ptr Q1, T1, Q2, T2;

        m = a + (b - a) / 2;

        n1 = m - a;
        n2 = b - m;
        q1len = FLINT_MIN(2 * n1 + 1, num);
        t1len = FLINT_MIN(2 * n1 - 1, num);
        q2len = FLINT_MIN(2 * n2 + 1, num);
        t2len = FLINT_MIN(2 * n2 - 1, num);
        qlen = FLINT_MIN(q1len + q2len - 1, num);
        tlen = FLINT_MIN(t1len + q2len - 1, num);

        alloc = q1len + q2len + t1len + t2len;
        Q1 = _arb_vec_init(alloc);
        Q2 = Q1 + q1len;
        T1 = Q2 + q2len;
        T2 = T1 + t1len;

        bsplit(Q1, T1, z, a, m, num, prec);
        bsplit(Q2, T2, z, m, b, num, prec);

        _arb_poly_mullow(Q, Q2, q2len, Q1, q1len, qlen, prec);
        _arb_poly_mullow(T, Q2, q2len, T1, t1len, tlen, prec);
        _arb_poly_add(T, T, tlen, T2, t2len, prec);

        _arb_vec_clear(Q1, alloc);
    }
}

void
_arb_poly_mullow_cpx(arb_ptr res, arb_srcptr src, slong len, const arb_t c, slong trunc, slong prec)
{
    slong i;

    if (len < trunc)
        arb_set(res + len, src + len - 1);

    for (i = len - 1; i > 0; i--)
    {
        arb_mul(res + i, src + i, c, prec);
        arb_add(res + i, res + i, src + i - 1, prec);
    }

    arb_mul(res, src, c, prec);
}

void
_arb_poly_log_cpx_series(arb_ptr res, const arb_t c, slong num, slong prec)
{
    slong i;

    for (i = 0; i < num; i++)
    {
        if (i == 0)
            arb_log(res + i, c, prec);
        else if (i == 1)
            arb_inv(res + i, c, prec);
        else
            arb_mul(res + i, res + i - 1, res + 1, prec);
    }

    for (i = 2; i < num; i++)
    {
        arb_div_ui(res + i, res + i, i, prec);

        if (i % 2 == 0)
            arb_neg(res + i, res + i);
    }
}

void
_arb_poly_gamma_stirling_eval2(arb_ptr res, const arb_t z, slong n, slong num, int diff, slong prec)
{
    slong k, tlen, qlen;
    arb_ptr T, Q;
    mag_ptr err;
    arb_t c;

    T = _arb_vec_init(num);
    Q = _arb_vec_init(num);
    err = _mag_vec_init(num);
    arb_init(c);

    arb_gamma_stirling_bound(err, z, 0, num, n);

    if (n <= 1)
    {
        _arb_vec_zero(res, num);
    }
    else
    {
        qlen = FLINT_MIN(2 * (n - 1) + 1, num);
        tlen = FLINT_MIN(2 * (n - 1) - 1, num);
        bsplit(Q, T, z, 1, n, num, prec);
        _arb_poly_div_series(res, T, tlen, Q, qlen, num, prec);
    }

    if (diff)
    {
        _arb_vec_add_error_mag_vec(res, err, num);
        _arb_poly_derivative(res, res, num, prec);

        if (num > 1)
        {
            /* add log(z+x) - 1/(2(z+x)) */
            arb_inv(c, z, prec);
            _arb_vec_set_powers(T, c, num, prec);

            for (k = 1; k < num - 1; k++)
            {
                arb_mul_2exp_si(T, z, 1);
                arb_div_ui(T, T, k, prec);
                arb_add_ui(T, T, 1, prec);
                arb_mul_2exp_si(T, T, -1);

                if (k % 2 == 0)
                    arb_submul(res + k, T, T + k + 1, prec);
                else
                    arb_addmul(res + k, T, T + k + 1, prec);
            }

            arb_mul_2exp_si(c, c, -1);
            arb_sub(res, res, c, prec);

            arb_log(c, z, prec);
            arb_add(res, res, c, prec);
        }
    }
    else
    {
        /* ((z-1/2) + t) * log(z+t) */
        _arb_poly_log_cpx_series(T, z, num, prec);
        arb_one(c);
        arb_mul_2exp_si(c, c, -1);
        arb_sub(c, z, c, prec);
        _arb_poly_mullow_cpx(T, T, num, c, num, prec);

        /* constant term */
        arb_const_log_sqrt2pi(c, prec);
        arb_add(T, T, c, prec);

        /* subtract (z+t) */
        arb_sub(T, T, z, prec);
        if (num > 1)
            arb_sub_ui(T + 1, T + 1, 1, prec);

        _arb_vec_add(res, res, T, num, prec);

        _arb_vec_add_error_mag_vec(res, err, num);
    }

    _arb_vec_clear(T, num);
    _arb_vec_clear(Q, num);
    _mag_vec_clear(err, num);
    arb_clear(c);
}

void
_arb_poly_gamma_stirling_eval(arb_ptr res, const arb_t z, slong n, slong num, slong prec)
{
    _arb_poly_gamma_stirling_eval2(res, z, n, num, 0, prec);
}

void
_arb_poly_gamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec)
{
    int reflect;
    slong i, rflen, r, n, wp;
    arb_ptr t, u, v;
    arb_struct f[2];

    if (hlen == 1)
    {
        arb_gamma(res, h, prec);
        if (arb_is_finite(res))
            _arb_vec_zero(res + 1, len - 1);
        else
            _arb_vec_indeterminate(res + 1, len - 1);
        return;
    }

    hlen = FLINT_MIN(hlen, len);
    wp = prec + FLINT_BIT_COUNT(prec);

    t = _arb_vec_init(len);
    u = _arb_vec_init(len);
    v = _arb_vec_init(len);
    arb_init(f);
    arb_init(f + 1);

    /* use zeta values at small integers */
    if (arb_is_int(h) && (arf_cmpabs_ui(arb_midref(h), prec / 2) < 0))
    {
        r = arf_get_si(arb_midref(h), ARF_RND_DOWN);

        if (r <= 0)
        {
            _arb_vec_indeterminate(v, len);
        }
        else if (r == 1)
        {
            _arb_poly_lgamma_series_at_one(u, len, wp);
            _arb_poly_exp_series(v, u, len, len, wp);
        }
        else
        {
            _arb_poly_lgamma_series_at_one(u, len, wp);
            _arb_poly_exp_series(t, u, len, len, wp);
            arb_one(f);
            arb_one(f + 1);
            rflen = FLINT_MIN(len, r);
            _arb_poly_rising_ui_series(u, f, FLINT_MIN(2, len), r - 1, rflen, wp);
            _arb_poly_mullow(v, t, len, u, rflen, len, wp);
        }
    }
    else
    {
        /* otherwise use Stirling series */
        arb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, h, 1, 0, wp);

        /* gamma(h) = (rf(1-h, r) * pi) / (gamma(1-h+r) sin(pi h)), h = h0 + t*/
        if (reflect)
        {
            /* u = 1/gamma(r+1-h) */
            arb_sub_ui(f, h, r + 1, wp);
            arb_neg(f, f);
            _arb_poly_gamma_stirling_eval(t, f, n, len, wp);
            _arb_vec_neg(t, t, len);
            _arb_poly_exp_series(u, t, len, len, wp);
            for (i = 1; i < len; i += 2)
                arb_neg(u + i, u + i);

            /* v = 1/sin(pi x) */
            arb_set(f, h);
            arb_one(f + 1);
            _arb_poly_sin_pi_series(t, f, 2, len, wp);
            _arb_poly_inv_series(v, t, len, len, wp);

            _arb_poly_mullow(t, u, len, v, len, len, wp);

            /* rf(1-h,r) * pi */
            if (r == 0)
            {
                rflen = 1;
                arb_const_pi(u, wp);
            }
            else
            {
                arb_sub_ui(f, h, 1, wp);
                arb_neg(f, f);
                arb_set_si(f + 1, -1);
                rflen = FLINT_MIN(len, r + 1);
                _arb_poly_rising_ui_series(u, f, FLINT_MIN(2, len), r, rflen, wp);
                arb_const_pi(v, wp);
                _arb_vec_scalar_mul(u, u, rflen, v, wp);
            }

            /* multiply by rising factorial */
            _arb_poly_mullow(v, t, len, u, rflen, len, wp);
        }
        else
        {
            /* gamma(h) = gamma(h+r) / rf(h,r) */
            if (r == 0)
            {
                arb_add_ui(f, h, r, wp);
                _arb_poly_gamma_stirling_eval(t, f, n, len, wp);
                _arb_poly_exp_series(v, t, len, len, wp);
            }
            else
            {
                /* TODO: div_series may be better (once it has a good basecase),
                         if the rising factorial is short */
                arb_set(f, h);
                arb_one(f + 1);
                rflen = FLINT_MIN(len, r + 1);
                _arb_poly_rising_ui_series(u, f, FLINT_MIN(2, len), r, rflen, wp);
                _arb_poly_inv_series(t, u, rflen, len, wp);

                arb_add_ui(f, h, r, wp);
                _arb_poly_gamma_stirling_eval(v, f, n, len, wp);
                _arb_poly_exp_series(u, v, len, len, wp);

                _arb_poly_mullow(v, u, len, t, len, len, wp);
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
arb_poly_gamma_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
{
    arb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _arb_vec_indeterminate(res->coeffs, n);
    else
        _arb_poly_gamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

