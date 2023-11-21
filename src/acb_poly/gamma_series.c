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

void acb_gamma_stirling_bound(mag_ptr err, const acb_t x, slong k0, slong knum, slong n);

void acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t x, int use_reflect, int digamma, slong prec);

void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);

static void
bsplit(acb_ptr Q, acb_ptr T, const acb_t z, slong a, slong b, slong num, slong prec)
{
    if (b - a == 1)
    {
        arb_gamma_stirling_coeff(acb_realref(T), a, 0, prec);
        arb_zero(acb_imagref(T));

        if (a == 1)
        {   /* (z + t) */
            acb_set(Q, z);
            if (num > 1) acb_one(Q + 1);
            if (num > 2) acb_zero(Q + 2);
        }
        else
        {   /* (z + t)^2 */
            acb_mul(Q, z, z, prec);  /* TODO: precompute */
            if (num > 1) acb_mul_2exp_si(Q + 1, z, 1);
            if (num > 2) acb_one(Q + 2);
        }
    }
    else
    {
        slong m, n1, n2, q1len, q2len, t1len, t2len, qlen, tlen, alloc;
        acb_ptr Q1, T1, Q2, T2;

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
        Q1 = _acb_vec_init(alloc);
        Q2 = Q1 + q1len;
        T1 = Q2 + q2len;
        T2 = T1 + t1len;

        bsplit(Q1, T1, z, a, m, num, prec);
        bsplit(Q2, T2, z, m, b, num, prec);

        _acb_poly_mullow(Q, Q2, q2len, Q1, q1len, qlen, prec);
        _acb_poly_mullow(T, Q2, q2len, T1, t1len, tlen, prec);
        _acb_poly_add(T, T, tlen, T2, t2len, prec);

        _acb_vec_clear(Q1, alloc);
    }
}

void
_acb_poly_log_cpx_series(acb_ptr res, const acb_t c, slong num, slong prec)
{
    slong i;

    for (i = 0; i < num; i++)
    {
        if (i == 0)
            acb_log(res + i, c, prec);
        else if (i == 1)
            acb_inv(res + i, c, prec);
        else
            acb_mul(res + i, res + i - 1, res + 1, prec);
    }

    for (i = 2; i < num; i++)
    {
        acb_div_ui(res + i, res + i, i, prec);

        if (i % 2 == 0)
            acb_neg(res + i, res + i);
    }
}

void
_acb_poly_gamma_stirling_eval2(acb_ptr res, const acb_t z, slong n, slong num, int diff, slong prec)
{
    slong k, tlen, qlen;
    acb_ptr T, Q;
    mag_ptr err;
    acb_t c;

    T = _acb_vec_init(num);
    Q = _acb_vec_init(num);
    err = _mag_vec_init(num);
    acb_init(c);

    acb_gamma_stirling_bound(err, z, 0, num, n);

    if (n <= 1)
    {
        _acb_vec_zero(res, num);
    }
    else
    {
        qlen = FLINT_MIN(2 * (n - 1) + 1, num);
        tlen = FLINT_MIN(2 * (n - 1) - 1, num);
        bsplit(Q, T, z, 1, n, num, prec);
        _acb_poly_div_series(res, T, tlen, Q, qlen, num, prec);
    }

    if (diff)
    {
        _acb_vec_add_error_mag_vec(res, err, num);
        _acb_poly_derivative(res, res, num, prec);

        if (num > 1)
        {
            /* add log(z+x) - 1/(2(z+x)) */
            acb_inv(c, z, prec);
            _acb_vec_set_powers(T, c, num, prec);

            for (k = 1; k < num - 1; k++)
            {
                acb_mul_2exp_si(T, z, 1);
                acb_div_ui(T, T, k, prec);
                acb_add_ui(T, T, 1, prec);
                acb_mul_2exp_si(T, T, -1);

                if (k % 2 == 0)
                    acb_submul(res + k, T, T + k + 1, prec);
                else
                    acb_addmul(res + k, T, T + k + 1, prec);
            }

            acb_mul_2exp_si(c, c, -1);
            acb_sub(res, res, c, prec);

            acb_log(c, z, prec);
            acb_add(res, res, c, prec);
        }
    }
    else
    {
        /* ((z-1/2) + t) * log(z+t) */
        _acb_poly_log_cpx_series(T, z, num, prec);
        acb_one(c);
        acb_mul_2exp_si(c, c, -1);
        acb_sub(c, z, c, prec);
        _acb_poly_mullow_cpx(T, T, num, c, num, prec);

        /* constant term */
        arb_const_log_sqrt2pi(acb_realref(c), prec);
        arb_zero(acb_imagref(c));
        acb_add(T, T, c, prec);

        /* subtract (z+t) */
        acb_sub(T, T, z, prec);
        if (num > 1)
            acb_sub_ui(T + 1, T + 1, 1, prec);

        _acb_vec_add(res, res, T, num, prec);

        _acb_vec_add_error_mag_vec(res, err, num);
    }

    _acb_vec_clear(T, num);
    _acb_vec_clear(Q, num);
    _mag_vec_clear(err, num);
    acb_clear(c);
}

void
_acb_poly_gamma_stirling_eval(acb_ptr res, const acb_t z, slong n, slong num, slong prec)
{
    _acb_poly_gamma_stirling_eval2(res, z, n, num, 0, prec);
}

void
_acb_poly_gamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong len, slong prec)
{
    int reflect;
    slong i, rflen, r, n, wp;
    acb_ptr t, u, v;
    acb_struct f[2];

    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_gamma(res, h, prec);
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
        _arb_poly_gamma_series(tmp, tmp, hlen, len, prec);
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

    /* use Stirling series */
    acb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, h, 1, 0, wp);

    /* gamma(h) = (rf(1-h, r) * pi) / (gamma(1-h+r) sin(pi h)), h = h0 + t*/
    if (reflect)
    {
        /* u = 1/gamma(r+1-h) */
        acb_sub_ui(f, h, r + 1, wp);
        acb_neg(f, f);
        _acb_poly_gamma_stirling_eval(t, f, n, len, wp);
        _acb_vec_neg(t, t, len);
        _acb_poly_exp_series(u, t, len, len, wp);
        for (i = 1; i < len; i += 2)
            acb_neg(u + i, u + i);

        /* v = 1/sin(pi x) */
        acb_set(f, h);
        acb_one(f + 1);
        _acb_poly_sin_pi_series(t, f, 2, len, wp);
        _acb_poly_inv_series(v, t, len, len, wp);

        _acb_poly_mullow(t, u, len, v, len, len, wp);

        /* rf(1-h,r) * pi */
        if (r == 0)
        {
            rflen = 1;
            acb_const_pi(u, wp);
        }
        else
        {
            acb_sub_ui(f, h, 1, wp);
            acb_neg(f, f);
            acb_set_si(f + 1, -1);
            rflen = FLINT_MIN(len, r + 1);
            _acb_poly_rising_ui_series(u, f, FLINT_MIN(2, len), r, rflen, wp);
            acb_const_pi(v, wp);
            _acb_vec_scalar_mul(u, u, rflen, v, wp);
        }

        /* multiply by rising factorial */
        _acb_poly_mullow(v, t, len, u, rflen, len, wp);
    }
    else
    {
        /* gamma(h) = gamma(h+r) / rf(h,r) */
        if (r == 0)
        {
            acb_add_ui(f, h, r, wp);
            _acb_poly_gamma_stirling_eval(t, f, n, len, wp);
            _acb_poly_exp_series(v, t, len, len, wp);
        }
        else
        {
            /* TODO: div_series may be better (once it has a good basecase),
                     if the rising factorial is short */
            acb_set(f, h);
            acb_one(f + 1);
            rflen = FLINT_MIN(len, r + 1);
            _acb_poly_rising_ui_series(u, f, FLINT_MIN(2, len), r, rflen, wp);
            _acb_poly_inv_series(t, u, rflen, len, wp);

            acb_add_ui(f, h, r, wp);
            _acb_poly_gamma_stirling_eval(v, f, n, len, wp);
            _acb_poly_exp_series(u, v, len, len, wp);

            _acb_poly_mullow(v, u, len, t, len, len, wp);
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
acb_poly_gamma_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
{
    acb_poly_fit_length(res, n);

    if (f->length == 0 || n == 0)
        _acb_vec_indeterminate(res->coeffs, n);
    else
        _acb_poly_gamma_series(res->coeffs, f->coeffs, f->length, n, prec);

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

