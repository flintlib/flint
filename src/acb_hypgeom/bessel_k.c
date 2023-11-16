/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_bessel_k_asymp(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    acb_t t, a, b, w;

    acb_init(t);
    acb_init(a);
    acb_init(b);
    acb_init(w);

    acb_one(a);
    acb_mul_2exp_si(a, a, -1);
    acb_add(a, a, nu, prec);

    acb_mul_2exp_si(b, nu, 1);
    acb_add_ui(b, b, 1, prec);

    acb_mul_2exp_si(w, z, 1);

    acb_hypgeom_u_asymp(t, a, b, w, -1, prec);

    if (!scaled)
    {
        acb_neg(a, z);
        acb_exp(a, a, prec);
        acb_mul(t, t, a, prec);
    }

    acb_mul_2exp_si(w, z, 1);
    acb_rsqrt(w, w, prec);
    acb_mul(res, t, w, prec);

    arb_const_sqrt_pi(acb_realref(w), prec);
    acb_mul_arb(res, res, acb_realref(w), prec);

    acb_clear(t);
    acb_clear(a);
    acb_clear(b);
    acb_clear(w);
}

void
acb_hypgeom_bessel_k_0f1_series(acb_poly_t res,
    const acb_poly_t nu, const acb_poly_t z,
    int scaled, slong len, slong prec)
{
    acb_poly_t s, u, A, B;
    acb_poly_struct b[2];
    arb_t c;
    slong wlen;
    int singular;

    acb_poly_init(s);
    acb_poly_init(u);
    acb_poly_init(A);
    acb_poly_init(B);
    acb_poly_init(b + 0);
    acb_poly_init(b + 1);
    arb_init(c);

    singular = (nu->length == 0) || acb_is_int(nu->coeffs);
    wlen = len + (singular != 0);

    /* A = (z/2)^nu, B = 1/A */
    acb_poly_scalar_mul_2exp_si(A, z, -1);
    acb_poly_pow_series(A, A, nu, wlen, prec);
    acb_poly_inv_series(B, A, wlen, prec);

    acb_poly_mullow(u, z, z, wlen, prec);
    acb_poly_scalar_mul_2exp_si(u, u, -2);

    acb_poly_one(b + 1);
    acb_poly_add_si(b + 0, nu, 1, prec);
    acb_hypgeom_pfq_series_direct(s, NULL, 0, b, 2, u, 1, -1, wlen, prec);
    acb_poly_mullow(A, A, s, wlen, prec);

    acb_poly_add_si(b + 0, nu, -1, prec);
    acb_poly_neg(b + 0, b + 0);
    acb_hypgeom_pfq_series_direct(s, NULL, 0, b, 2, u, 1, -1, wlen, prec);
    acb_poly_mullow(B, B, s, wlen, prec);

    acb_poly_sub(A, B, A, prec);
    acb_poly_scalar_mul_2exp_si(A, A, -1);

    /* multiply by pi csc(pi nu) */
    acb_poly_sin_pi_series(B, nu, wlen, prec);

    if (singular)
    {
        acb_poly_shift_right(A, A, 1);
        acb_poly_shift_right(B, B, 1);
    }

    if (scaled)
    {
        acb_poly_exp_series(u, z, len, prec);
        acb_poly_mullow(A, A, u, len, prec);
    }

    acb_poly_div_series(res, A, B, len, prec);

    arb_const_pi(c, prec);
    _acb_vec_scalar_mul_arb(res->coeffs, res->coeffs, res->length, c, prec);

    acb_poly_clear(s);
    acb_poly_clear(u);
    acb_poly_clear(A);
    acb_poly_clear(B);
    acb_poly_clear(b + 0);
    acb_poly_clear(b + 1);
    arb_clear(c);
}

void
acb_hypgeom_bessel_k_0f1(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    if (acb_is_int(nu))
    {
        acb_poly_t nux, zx, rx;

        acb_poly_init(nux);
        acb_poly_init(zx);
        acb_poly_init(rx);

        acb_poly_set_coeff_acb(nux, 0, nu);
        acb_poly_set_coeff_si(nux, 1, 1);
        acb_poly_set_acb(zx, z);

        acb_hypgeom_bessel_k_0f1_series(rx, nux, zx, scaled, 1, prec);

        acb_poly_get_coeff_acb(res, rx, 0);

        acb_poly_clear(nux);
        acb_poly_clear(zx);
        acb_poly_clear(rx);
    }
    else
    {
        acb_t t, u, v, w;
        acb_struct b[2];

        acb_init(t);
        acb_init(u);
        acb_init(v);
        acb_init(w);
        acb_init(b + 0);
        acb_init(b + 1);

        /* u = 0F1(1+nu), v = 0F1(1-nu) */
        acb_mul(t, z, z, prec);
        acb_mul_2exp_si(t, t, -2);
        acb_add_ui(b, nu, 1, prec);
        acb_one(b + 1);
        acb_hypgeom_pfq_direct(u, NULL, 0, b, 2, t, -1, prec);
        acb_sub_ui(b, nu, 1, prec);
        acb_neg(b, b);
        acb_hypgeom_pfq_direct(v, NULL, 0, b, 2, t, -1, prec);

        /* v = v * gamma(nu) / (z/2)^nu */
        acb_mul_2exp_si(t, z, -1);
        acb_pow(t, t, nu, prec);
        acb_gamma(w, nu, prec);
        acb_mul(v, v, w, prec);
        acb_div(v, v, t, prec);

        /* u = u * t * pi / (gamma(nu) * nu * sin(pi nu)) */
        acb_mul(u, u, t, prec);
        acb_const_pi(t, prec);
        acb_mul(u, u, t, prec);
        acb_sin_pi(t, nu, prec);
        acb_mul(t, t, w, prec);
        acb_mul(t, t, nu, prec);
        acb_div(u, u, t, prec);

        acb_sub(v, v, u, prec);
        acb_mul_2exp_si(v, v, -1);

        if (scaled)
        {
            acb_exp(t, z, prec);
            acb_mul(v, v, t, prec);
        }

        acb_set(res, v);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
        acb_clear(w);
        acb_clear(b + 0);
        acb_clear(b + 1);
    }
}

void
acb_hypgeom_bessel_k_nointegration(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    mag_t zmag;

    mag_init(zmag);
    acb_get_mag(zmag, z);

    if (mag_cmp_2exp_si(zmag, 4) < 0 ||
        (mag_cmp_2exp_si(zmag, 64) < 0 && 2 * mag_get_d(zmag) < prec))
        acb_hypgeom_bessel_k_0f1(res, nu, z, scaled, prec);
    else
        acb_hypgeom_bessel_k_asymp(res, nu, z, scaled, prec);

    mag_clear(zmag);
}

void
_acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
{
    acb_t res2;
    slong acc, max, t;

    acb_init(res2);

    acb_hypgeom_bessel_k_nointegration(res2, nu, z, scaled, prec);

    acc = acb_rel_accuracy_bits(res2);

    if (acc < 0.5 * prec)
    {
        max = prec;
        t = acb_rel_accuracy_bits(z);
        max = FLINT_MIN(max, t);
        t = acb_rel_accuracy_bits(nu);
        max = FLINT_MIN(max, t);

        if (max > 2 && acc < 0.5 * max)
        {
            if (acb_is_real(nu) && acb_is_real(z) && arf_cmp_d(arb_midref(acb_realref(nu)), -0.5) > 0 &&
                arf_cmp_2exp_si(arb_midref(acb_realref(z)), -16) > 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 60) < 0)
            {
                arb_hypgeom_bessel_k_integration(acb_realref(res),
                        acb_realref(nu), acb_realref(z), scaled, prec);
                arb_zero(acb_imagref(res));

                if (acb_rel_accuracy_bits(res) > acb_rel_accuracy_bits(res2) ||
                    (acb_is_finite(res) && !acb_is_finite(res2)))
                {
                    acb_swap(res, res2);
                }
            }
        }
    }

    acb_swap(res, res2);
    acb_clear(res2);
}

void
acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    _acb_hypgeom_bessel_k(res, nu, z, 0, prec);
}

void
acb_hypgeom_bessel_k_scaled(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    _acb_hypgeom_bessel_k(res, nu, z, 1, prec);
}
