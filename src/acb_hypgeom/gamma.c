/*
    Copyright (C) 2014, 2015, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"
#include "arb_hypgeom.h"

void
acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t z, int use_reflect, int digamma, slong prec);

void acb_gamma_stirling_bound(mag_ptr err, const acb_t x, slong k0, slong knum, slong n);

void
acb_hypgeom_gamma_stirling_inner(acb_t s, const acb_t z, slong N, slong prec)
{
    acb_t logz, t;
    mag_t err;

    mag_init(err);
    acb_init(t);
    acb_init(logz);

    acb_gamma_stirling_bound(err, z, 0, 1, N);

    /* t = (z-0.5)*log(z) - z + log(2*pi)/2 */
    acb_log(logz, z, prec);
    arb_one(acb_realref(t));
    arb_mul_2exp_si(acb_realref(t), acb_realref(t), -1);
    acb_sub(t, z, t, prec);
    acb_mul(t, logz, t, prec);
    acb_sub(t, t, z, prec);
    arb_const_log_sqrt2pi(acb_realref(logz), prec);
    arb_add(acb_realref(t), acb_realref(t), acb_realref(logz), prec);

    /* sum part */
    if (prec <= 128 || (prec <= 1024 && N <= 40) || (prec <= 2048 && N <= 16))
        acb_hypgeom_gamma_stirling_sum_horner(s, z, N, prec);
    else
        acb_hypgeom_gamma_stirling_sum_improved(s, z, N, 0, prec);

    acb_add(s, s, t, prec);
    acb_add_error_mag(s, err);

    acb_clear(t);
    acb_clear(logz);
    mag_clear(err);
}


void
acb_hypgeom_gamma_stirling(acb_t y, const acb_t x, int reciprocal, slong prec)
{
    int reflect;
    slong r, n, wp;
    acb_t t, u, v;
    double acc;

    wp = prec + FLINT_BIT_COUNT(prec);

    /* todo: for large x (if exact or accurate enough), increase precision */
    acc = acb_rel_accuracy_bits(x);
    acc = FLINT_MAX(acc, 0);
    wp = FLINT_MIN(prec, acc + 20);
    wp = FLINT_MAX(wp, 2);
    wp = wp + FLINT_BIT_COUNT(wp);

    if (acc < 3)  /* try to avoid divisions blowing up */
    {
        if (arf_cmp_d(arb_midref(acb_realref(x)), -0.5) < 0)
        {
            reflect = 1;
            r = 0;
        }
        else if (arf_cmp_si(arb_midref(acb_realref(x)), 1) < 0)
        {
            reflect = 0;
            r = 1;
        }
        else
        {
            reflect = 0;
            r = 0;
        }

        n = 1;
    }
    else
    {
        acb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);
    }

    acb_init(t);
    acb_init(u);
    acb_init(v);

    if (reflect)
    {
        acb_sub_ui(t, x, 1, wp);
        acb_neg(t, t);
        acb_hypgeom_rising_ui_rec(u, t, r, wp);
        arb_const_pi(acb_realref(v), wp);
        acb_mul_arb(u, u, acb_realref(v), wp);
        acb_add_ui(t, t, r, wp);
        acb_hypgeom_gamma_stirling_inner(v, t, n, wp);

        if (reciprocal)
        {
            /* rgamma(x) = gamma(1-x+r) sin(pi x) / ((rf(1-x, r) * pi) */
            acb_exp(v, v, wp);
            acb_sin_pi(t, x, wp);
            acb_mul(v, v, t, wp);
            acb_mul(y, u, v, wp);
            acb_div(y, v, u, prec);
        }
        else
        {
            /* gamma(x) = (rf(1-x, r) * pi) rgamma(1-x+r) csc(pi x) */
            acb_neg(v, v);
            acb_exp(v, v, wp);
            acb_csc_pi(t, x, wp);
            acb_mul(v, v, t, wp);
            acb_mul(y, v, u, prec);
        }
    }
    else
    {
        acb_add_ui(t, x, r, wp);
        acb_hypgeom_gamma_stirling_inner(u, t, n, wp);

        if (reciprocal)
        {
            /* rgamma(x) = rf(x,r) rgamma(x+r) */
            acb_neg(u, u);
            acb_exp(u, u, prec);
            acb_hypgeom_rising_ui_rec(v, x, r, wp);
            acb_mul(y, v, u, prec);
        }
        else
        {
            /* gamma(x) = gamma(x+r) / rf(x,r) */
            acb_exp(u, u, prec);
            acb_hypgeom_rising_ui_rec(v, x, r, wp);
            acb_div(y, u, v, prec);
        }
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_hypgeom_gamma(acb_t y, const acb_t x, slong prec)
{
    if (acb_is_real(x))
    {
        arb_hypgeom_gamma(acb_realref(y), acb_realref(x), prec);
        arb_zero(acb_imagref(y));
        return;
    }

    if (acb_hypgeom_gamma_taylor(y, x, 0, prec))
        return;

    acb_hypgeom_gamma_stirling(y, x, 0, prec);
}

void
acb_hypgeom_rgamma(acb_t y, const acb_t x, slong prec)
{
    mag_t magz;

    if (acb_is_real(x))
    {
        arb_hypgeom_rgamma(acb_realref(y), acb_realref(x), prec);
        arb_zero(acb_imagref(y));
        return;
    }

    if (acb_hypgeom_gamma_taylor(y, x, 1, prec))
        return;

    mag_init(magz);
    acb_get_mag(magz, x);

    if (mag_is_inf(magz))
    {
        acb_indeterminate(y);
    }
    else
    {
        acb_hypgeom_gamma_stirling(y, x, 1, prec);

        /* Todo: improved bounds computation */
        if (!acb_is_finite(y))
        {
            arb_t t, u, R;

            arb_init(R);
            arb_init(t);
            arb_init(u);

            arf_set_mag(arb_midref(R), magz);

            arb_set_d(u, 0.5);
            arb_add(u, u, R, MAG_BITS);
            arb_pow(u, R, u, MAG_BITS);

            arb_const_pi(t, MAG_BITS);
            arb_mul(t, t, R, MAG_BITS);
            arb_mul_2exp_si(t, t, -1);
            arb_exp(t, t, MAG_BITS);

            arb_mul(t, t, u, MAG_BITS);

            arb_get_mag(magz, t);

            acb_zero(y);
            acb_add_error_mag(y, magz);

            arb_clear(R);
            arb_clear(t);
            arb_clear(u);
        }
    }

    mag_clear(magz);
}

