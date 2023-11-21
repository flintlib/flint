/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

/* note: will not return a wrong value, as arf_get_si aborts on overflow */
slong
arb_get_si_lower(const arb_t x)
{
    arf_t t;
    slong v;

    arf_init(t);
    arf_set_mag(t, arb_radref(x));
    arf_sub(t, arb_midref(x), t, 2 * FLINT_BITS, ARF_RND_FLOOR);

    v = arf_get_si(t, ARF_RND_FLOOR);

    arf_clear(t);

    return v;
}

slong
polylog_choose_terms(mag_t err, slong sigma, const mag_t z, slong d, slong prec)
{
    slong N;

    for (N = 3; ; N = FLINT_MAX(N+3, N*1.1))
    {
        mag_polylog_tail(err, z, sigma, d, N);

        /* TODO: do something else when |Li_s(z)| is very small/very large? */
        if (mag_cmp_2exp_si(err, -prec) < 0)
            break;

        if (N > 100 * prec)
        {
            N = 3;
            mag_inf(err);
            break;
        }
    }

    return N;
}

int
polylog_is_real(const acb_t s, const acb_t z)
{
    if (!arb_is_zero(acb_imagref(s)))
        return 0;
    else if (!arb_is_zero(acb_imagref(z)))
        return 0;
    else if (arb_contains_si(acb_realref(z), 1))
        return 0;
    else if (acb_is_int(s) && arb_is_nonpositive(acb_realref(s)))
        return 1;
    else
        return (arf_cmp_2exp_si(arb_midref(acb_realref(z)), 0) < 0);
}

void
_acb_poly_polylog_cpx_zeta(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec)
{
    acb_ptr e1, e2, z1, z2, e1z1, e2z2;
    acb_t t, u, v;
    slong k, len2;
    int deflate_zeta, deflate_gamma, is_real;

    if (!acb_is_finite(s) || !acb_is_finite(z))
    {
        _acb_vec_indeterminate(w, len);
        return;
    }

    if (acb_is_one(z))
    {
        if (arb_gt(acb_realref(s), acb_realref(z))) /* Re(s) > 1 */
        {
            acb_zeta(w, s, prec);
            _acb_vec_indeterminate(w + 1, len - 1);
        }
        else
        {
            _acb_vec_indeterminate(w, len);
        }

        return;
    }

    is_real = polylog_is_real(s, z);

    acb_init(t);
    acb_init(u);
    acb_init(v);

    /* v = 1-s */
    acb_one(v);
    acb_sub(v, v, s, prec);

    /* pole of zeta */
    deflate_zeta = acb_is_one(v);

    /* poles of gamma at nonpositive integer v */
    deflate_gamma = (arb_is_zero(acb_imagref(v)) &&
            arb_is_int(acb_realref(v)) &&
            arf_sgn(arb_midref(acb_realref(v))) <= 0);

    len2 = len + deflate_gamma;

    e1 = _acb_vec_init(len + 1);
    e2 = _acb_vec_init(len + 1);
    z1 = _acb_vec_init(len + 1);
    z2 = _acb_vec_init(len + 1);
    e1z1 = _acb_vec_init(len + 1);
    e2z2 = _acb_vec_init(len + 1);

    /* u = log(-z)/(pi*i) */
    acb_neg(t, z);
    acb_log(t, t, prec);
    acb_const_pi(u, prec);
    acb_mul_onei(u, u);
    acb_div(u, t, u, prec);

    /* z1 = zeta(v+x, 1/2 + log(-z)/(2*pi*i)) */
    acb_one(t);
    acb_add(t, t, u, prec);
    acb_mul_2exp_si(t, t, -1);
    _acb_poly_zeta_cpx_series(z1, v, t, deflate_zeta, len2, prec);

    /* z2 = zeta(v+x, 1/2 - log(-z)/(2*pi*i)) */
    acb_one(t);
    acb_sub(t, t, u, prec);
    acb_mul_2exp_si(t, t, -1);
    _acb_poly_zeta_cpx_series(z2, v, t, deflate_zeta, len2, prec);

    /* e1 = (i/(2pi))^(v+x) */
    acb_onei(t);
    acb_const_pi(u, prec);
    acb_div(t, t, u, prec);
    acb_mul_2exp_si(t, t, -1);
    _acb_poly_acb_pow_cpx(e1, t, v, len + (deflate_zeta || deflate_gamma), prec);

    /* e2 = (1/(2 pi i))^(v+x) */
    acb_conj(t, t);
    _acb_poly_acb_pow_cpx(e2, t, v, len + (deflate_zeta || deflate_gamma), prec);

    _acb_poly_mullow(e1z1, e1, len2, z1, len2, len2, prec);
    _acb_poly_mullow(e2z2, e2, len2, z2, len2, len2, prec);
    _acb_vec_add(z1, e1z1, e2z2, len2, prec);

    if (deflate_gamma)
    {
        /* gamma(v+x) = pi/sin(pi(v+x)) * 1/gamma(1-v-x) */

        /* TODO: write a csc function? */
        acb_zero(e1);
        acb_const_pi(e1 + 1, prec);
        acb_mul_2exp_si(e2, v, -1);
        if (!arb_is_int(acb_realref(e2)))
            acb_neg(e1 + 1, e1 + 1);
        _acb_poly_sin_series(e2, e1, 2, len2, prec);
        _acb_poly_inv_series(e1, e2 + 1, len, len, prec);
        acb_const_pi(e2, prec);
        _acb_vec_scalar_mul(e1, e1, len, e2, prec);

        acb_set(z2, s);
        acb_set_si(z2 + 1, -1);
        _acb_poly_rgamma_series(e2, z2, 2, len, prec);
        _acb_poly_mullow(z2, e1, len, e2, len, len, prec);

        _acb_poly_mullow(w, z1 + 1, len, z2, len, len, prec);
    }
    else
    {
        if (deflate_zeta)
        {
            for (k = 0; k < len; k++)
            {
                arb_mul_2exp_si(acb_realref(e1 + k + 1), acb_realref(e1 + k + 1), 1);
                arb_add(acb_realref(z1 + k), acb_realref(z1 + k), acb_realref(e1 + k + 1), prec);
            }

        }

        /* gamma(v+x) */
        acb_set(e1, v);
        if (len > 1)
            acb_one(e1 + 1);
        _acb_poly_gamma_series(z2, e1, FLINT_MIN(len, 2), len, prec);

        _acb_poly_mullow(w, z2, len, z1, len, len, prec);
    }

    /* correct signs (from s -> 1-s) */
    for (k = 1; k < len; k += 2)
        acb_neg(w + k, w + k);

    if (is_real)
        if (acb_is_finite(w))
            arb_zero(acb_imagref(w));

    _acb_vec_clear(e1, len + 1);
    _acb_vec_clear(e2, len + 1);
    _acb_vec_clear(z1, len + 1);
    _acb_vec_clear(z2, len + 1);
    _acb_vec_clear(e1z1, len + 1);
    _acb_vec_clear(e2z2, len + 1);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
_acb_poly_polylog_cpx_small(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec)
{
    slong k, N, sigma;
    int is_real;
    mag_t zmag, err, errf;
    acb_t a;

    acb_init(a);
    mag_init(zmag);
    mag_init(err);
    mag_init(errf);

    is_real = polylog_is_real(s, z);
    acb_get_mag(zmag, z);
    sigma = arb_get_si_lower(acb_realref(s));

    N = polylog_choose_terms(err, sigma, zmag, len - 1, prec);

    /* TODO: allow threading */
    acb_one(a);
    _acb_poly_powsum_series_naive(w, s, a, z, N - 1, len, prec);
    _acb_vec_scalar_mul(w, w, len, z, prec);

    for (k = 0; k < len; k++)
    {
        mag_polylog_tail(err, zmag, sigma, k, N);
        mag_rfac_ui(errf, k);
        mag_mul(err, err, errf);

        if (is_real && mag_is_finite(err))
            arb_add_error_mag(acb_realref(w + k), err);
        else
            acb_add_error_mag(w + k, err);
    }

    acb_clear(a);
    mag_clear(zmag);
    mag_clear(err);
    mag_clear(errf);
}

void
_acb_poly_polylog_cpx(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec)
{
    mag_t zmag;

    if (len == 1 && acb_equal_si(s, 2))
    {
        acb_hypgeom_dilog(w, z, prec);
        return;
    }

    mag_init(zmag);
    acb_get_mag(zmag, z);

    if (mag_cmp_2exp_si(zmag, -1) < 0)
        _acb_poly_polylog_cpx_small(w, s, z, len, prec);
    else
        _acb_poly_polylog_cpx_zeta(w, s, z, len, prec);

    mag_clear(zmag);
}

void
_acb_poly_polylog_series(acb_ptr res, acb_srcptr s, slong slen, const acb_t z, slong len, slong prec)
{
    acb_ptr t, u;

    slen = FLINT_MIN(slen, len);

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);

    _acb_poly_polylog_cpx(t, s, z, len, prec);

    /* compose with nonconstant part */
    acb_zero(u);
    _acb_vec_set(u + 1, s + 1, slen - 1);
    _acb_poly_compose_series(res, t, len, u, slen, len, prec);

    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
}

void
acb_poly_polylog_series(acb_poly_t res, const acb_poly_t s, const acb_t z, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, n);

    if (s->length == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_poly_polylog_series(res->coeffs, t, 1, z, n, prec);
        acb_clear(t);
    }
    else
    {
        _acb_poly_polylog_series(res->coeffs, s->coeffs, s->length, z, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}
