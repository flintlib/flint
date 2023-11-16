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

#ifdef __GNUC__
# define fabs __builtin_fabs
# define floor __builtin_floor
#else
# include <math.h>
#endif

void acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t z, int use_reflect, int digamma, slong prec);

void acb_hypgeom_gamma_stirling_inner(acb_t s, const acb_t z, slong N, slong prec);

static double
want_taylor(double x, double y, slong prec)
{
    if (y < 0.0) y = -y;
    if (x < 0.0) x = -2.0 * x;

    if ((prec < 128 && y > 4.0) || (prec < 256 && y > 5.0) ||
        (prec < 512 && y > 8.0) || (prec < 1024 && y > 9.0) || y > 10.0)
    {
        return 0;
    }

    if (x * (1.0 + 0.75 * y) > 8 + 0.15 * prec)
    {
        return 0;
    }

    return 1;
}

/* Linear fit on [0.5, 1.5] for
    lambda x: findroot(lambda y: im(loggamma(x+1j*y)) - (n+0.5)*pi */

static const double Atab[] = {
    4.5835631239879990091,
    6.4037921417161376741,
    7.9938623618272375768,
    9.4449131928216797873,
    10.802608819487725856,
    12.0918817314347272,
};

static const double Btab[] = {
    -1.1432582881376479127,
    -0.86248117216701645437,
    -0.75778990135448922722,
    -0.69734688055939976228,
    -0.65626499937495627271,
    -0.62578331900739100617,
};

void
_arb_const_log_pi(arb_t t, slong prec)
{
    arb_const_pi(t, prec + 2);
    arb_log(t, t, prec);
}

ARB_DEF_CACHED_CONSTANT(arb_const_log_pi, _arb_const_log_pi)

int
acb_hypgeom_lgamma_taylor(acb_t res, const acb_t z, slong prec)
{
    double x, y, acc;
    slong k, r, wp;
    acb_t t, u;
    int reflect;

    /* Assume xerr, yerr <= 1/16 */
    if (mag_cmp_2exp_si(arb_radref(acb_realref(z)), -4) > 0)
        return 0;

    if (mag_cmp_2exp_si(arb_radref(acb_imagref(z)), -4) > 0)
        return 0;

    acc = acb_rel_accuracy_bits(z);
    acc = FLINT_MAX(acc, 0);
    wp = FLINT_MIN(prec, acc + 20);
    wp = FLINT_MAX(wp, 2);

    /* x, y plus eventual rounding error */
    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);

    if (!want_taylor(x, y, wp))
        return 0;

    acb_init(t);
    acb_init(u);

    /* Reduce real part to (approximately) [0.5, 1.5]. */
    r = floor(x - 0.5);

    /* Reflection formula is slower but improves accuracy. */
    reflect = (x < -3.0);

    if (reflect)
    {
        acb_neg(u, z);
        acb_add_si(u, u, 2 + r, 2 * prec + 10);
        x = 2.0 + r - x;
        y = -y;
    }
    else
    {
        acb_sub_si(u, z, r, 2 * prec + 10);
        x = x - r;
    }

    for (k = 0; k < 6; k++)
    {
        if (fabs(y) <= Atab[k] + Btab[k] * x)
        {
            if (!acb_hypgeom_gamma_taylor(t, u, 1, wp))
            {
                acb_clear(t);
                acb_clear(u);
                return 0;
            }

            if (k % 2 == 0)
            {
                acb_log(t, t, wp);
                acb_neg(t, t);
            }
            else
            {
                acb_neg(t, t);
                acb_log(t, t, wp);
                acb_neg(t, t);
            }

            if (k != 0)
            {
                arb_t pi;
                arb_init(pi);
                arb_const_pi(pi, wp);
                arb_addmul_si(acb_imagref(t), pi, (y > 0) ? k : -k, wp);
                arb_clear(pi);
            }

            if (reflect)
            {
                acb_t v;
                acb_init(v);

                /* loggamma(x) = log(pi) - lsin(x) - loggamma(2+r-x) - logrf(2+r-x, -r-1) */

                acb_hypgeom_log_rising_ui(v, u, -r-1, wp);
                acb_log_sin_pi(res, z, wp);
                acb_add(res, res, v, wp);
                acb_add(res, res, t, wp);
                acb_neg(res, res);

                arb_const_log_pi(acb_realref(t), wp);
                arb_zero(acb_imagref(t));
                acb_add(res, res, t, prec);

                acb_clear(v);
            }
            else if (r == 0)
            {
                acb_set_round(res, t, prec);
            }
            else if (r > 0)
            {
                acb_hypgeom_log_rising_ui(res, u, r, wp);
                acb_add(res, res, t, prec);
            }
            else
            {
                acb_hypgeom_log_rising_ui(res, z, -r, wp);
                acb_sub(res, t, res, prec);
            }

            acb_clear(t);
            acb_clear(u);
            return 1;
        }
    }

    acb_clear(t);
    acb_clear(u);
    return 0;
}

void
acb_hypgeom_lgamma(acb_t y, const acb_t x, slong prec)
{
    int reflect;
    slong r, n, wp;
    acb_t t, u, v;
    double acc;

    if (acb_is_real(x) && arb_is_positive(acb_realref(x)))
    {
        arb_hypgeom_lgamma(acb_realref(y), acb_realref(x), prec);
        arb_zero(acb_imagref(y));
        return;
    }

    if (acb_hypgeom_lgamma_taylor(y, x, prec))
        return;

    acc = acb_rel_accuracy_bits(x);
    acc = FLINT_MAX(acc, 0);
    wp = FLINT_MIN(prec, acc + 20);
    wp = FLINT_MAX(wp, 2);
    wp = wp + FLINT_BIT_COUNT(wp);

    acb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);

    acb_init(t);
    acb_init(u);
    acb_init(v);

    if (reflect)
    {
        /* log gamma(x) = log rf(1-x, r) - log gamma(1-x+r) - log sin(pi x) + log(pi) */
        acb_sub_ui(u, x, 1, wp);
        acb_neg(u, u);

        acb_hypgeom_log_rising_ui(t, u, r, wp);

        acb_add_ui(u, u, r, wp);
        acb_hypgeom_gamma_stirling_inner(v, u, n, wp);
        acb_sub(t, t, v, wp);

        acb_log_sin_pi(u, x, wp);
        acb_sub(t, t, u, wp);

        arb_const_log_pi(acb_realref(u), wp);
        arb_zero(acb_imagref(u));

        acb_add(y, t, u, wp);
    }
    else
    {
        /* log gamma(x) = log gamma(x+r) - log rf(x,r) */
        acb_add_ui(t, x, r, wp);
        acb_hypgeom_gamma_stirling_inner(u, t, n, wp);
        acb_hypgeom_log_rising_ui(t, x, r, wp);
        acb_sub(y, u, t, prec);
    }

    if (!acb_is_finite(y))
        acb_indeterminate(y);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

