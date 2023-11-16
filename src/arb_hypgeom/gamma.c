/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "bernoulli.h"

/* tuning factor */
double GAMMA_STIRLING_BETA = 0.0;

#define PI 3.1415926535897932385

static slong
choose_n(double log2z, double argz, int digamma, slong prec)
{
    double argf, boundn, boundn_best;
    slong n, nbest;

    argf = 1.0 / cos(0.5 * argz);
    argf = log(argf) * (1. / log(2));

    boundn_best = 1e300;
    nbest = 1;

    for (n = 1; ; n++)
    {
        if (digamma)
            boundn = bernoulli_bound_2exp_si(2*n) - (2*n)*log2z + (2*n+1)*argf;
        else
            boundn = bernoulli_bound_2exp_si(2*n) - (2*n-1)*log2z + (2*n)*argf;

        /* success */
        if (boundn <= -prec)
            return n;

        if (boundn < boundn_best)
        {
            nbest = n;
            boundn_best = boundn;
        }

        /* if the term magnitude does not decrease, r is too small */
        if (boundn > 1)
        {
            /* printf("failure: prec = %ld, nbound_best = %f [%ld, %ld]\n", prec, boundn_best, n, nbest); */
            return nbest;
        }
    }
}

static void
choose_small(int * reflect, slong * r, slong * n,
    double x, double y, int use_reflect, int digamma, slong prec)
{
    double w, argz, log2z, BETA;
    slong rr;

    /* use reflection formula if very negative */
    if (x < -5.0 && use_reflect)
    {
        *reflect = 1;
        x = 1.0 - x;
    }
    else
    {
        *reflect = 0;
    }

    BETA = GAMMA_STIRLING_BETA;

    if (BETA < 0.12)
    {
        if (prec <= 32768)
            BETA = 0.17;
        else if (prec <= 131072)
            BETA = 0.20;
        else
            BETA = 0.24;
    }

    /* argument reduction until |z| >= w */
    w = FLINT_MAX(1.0, BETA * prec);

    rr = 0;
    while (x < 1.0 || x*x + y*y < w*w)
    {
        x++;
        rr++;
    }

    log2z = 0.5 * log(x*x + y*y) * 1.44269504088896341;
    argz = atan2(y, x);

    *r = rr;
    *n = choose_n(log2z, argz, digamma, prec);
}

static void
choose_large(int * reflect, slong * r, slong * n,
    const arf_t a, const arf_t b, int use_reflect, int digamma, slong prec)
{
    if (use_reflect && arf_sgn(a) < 0)
        *reflect = 1;
    else
        *reflect = 0;

    *r = 0;

    /* so big that we will certainly have n = 0 */
    if (arf_cmpabs_2exp_si(a, WORD_MAX / 8) >= 0 ||
        arf_cmpabs_2exp_si(b, WORD_MAX / 8) >= 0)
    {
        *n = 0;
    }
    else
    {
        slong ab, bb;
        double log2z, argz;

        ab = arf_abs_bound_lt_2exp_si(a);
        bb = arf_abs_bound_lt_2exp_si(b);

        log2z = FLINT_MAX(ab, bb);

        /* piecewise approximation of the argument */
        if (arf_is_zero(b))
        {
            if ((arf_sgn(a) < 0) && !(*reflect))
                argz = PI;
            else
                argz = 0.0;
        }
        else
        {
            if ((arf_sgn(a) < 0) && !(*reflect))
                if (arf_cmpabs(a, b) <= 0)
                    argz = PI * 0.75;
                else
                    argz = PI;
            else
                if (arf_cmpabs(a, b) <= 0)
                    argz = PI * 0.25;
                else
                    argz = PI * 0.5;
        }

        if (argz == PI)
            *n = 0;
        else
            *n = choose_n(log2z, argz, digamma, prec);
    }
}


void
acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const acb_t z, int use_reflect, int digamma, slong prec)
{
    const arf_struct * a = arb_midref(acb_realref(z));
    const arf_struct * b = arb_midref(acb_imagref(z));

    if (!arf_is_finite(a) || !arf_is_finite(b))
    {
        *reflect = *r = *n = 0;
    }
    else if (arf_cmpabs_2exp_si(a, 40) > 0 || arf_cmpabs_2exp_si(b, 40) > 0)
    {
        choose_large(reflect, r, n, a, b, use_reflect, digamma, prec);
    }
    else
    {
        choose_small(reflect, r, n,
            arf_get_d(a, ARF_RND_UP),
            arf_get_d(b, ARF_RND_UP), use_reflect, digamma, prec);
    }
}

void
arb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n,
    const arb_t x, int use_reflect, int digamma, slong prec)
{
    const arf_struct * a = arb_midref(x);

    if (arf_is_inf(a) || arf_is_nan(a))
    {
        *reflect = *r = *n = 0;
    }
    else if (arf_cmpabs_2exp_si(a, 40) > 0)
    {
        arf_t b;
        arf_init(b);
        choose_large(reflect, r, n, a, b, use_reflect, digamma, prec);
        arf_clear(b);
    }
    else
    {
        choose_small(reflect, r, n,
            arf_get_d(a, ARF_RND_UP), 0.0, use_reflect, digamma, prec);
    }
}

void arb_gamma_stirling_bound(mag_ptr err, const arb_t x, slong k0, slong knum, slong n);

void
arb_hypgeom_gamma_stirling_inner(arb_t s, const arb_t z, slong N, slong prec)
{
    arb_t logz, t;
    mag_t err;

    mag_init(err);
    arb_init(t);
    arb_init(logz);

    arb_gamma_stirling_bound(err, z, 0, 1, N);

    /* t = (z-0.5)*log(z) - z + log(2*pi)/2 */
    arb_log(logz, z, prec);
    arb_one(t);
    arb_mul_2exp_si(t, t, -1);
    arb_sub(t, z, t, prec);
    arb_mul(t, logz, t, prec);
    arb_sub(t, t, z, prec);
    arb_const_log_sqrt2pi(logz, prec);
    arb_add(t, t, logz, prec);

    /* sum part */
    if (prec <= 128 || (prec <= 768 && N <= 40) || (prec <= 2048 && N <= 16))
        arb_hypgeom_gamma_stirling_sum_horner(s, z, N, prec);
    else
        arb_hypgeom_gamma_stirling_sum_improved(s, z, N, 0, prec);

    arb_add(s, s, t, prec);

    mag_add(arb_radref(s), arb_radref(s), err);

    arb_clear(t);
    arb_clear(logz);
    mag_clear(err);
}

int
arb_hypgeom_gamma_exact(arb_t res, const arb_t x, int reciprocal, slong prec)
{
    if (arb_is_exact(x))
    {
        const arf_struct * mid = arb_midref(x);

        if (arf_is_special(mid))
        {
            if (!reciprocal && arf_is_pos_inf(mid))
                arb_set(res, x);
            else if (arf_is_nan(mid) || arf_is_neg_inf(mid) || !reciprocal)
                arb_indeterminate(res);
            else
                arb_zero(res);
            return 1;
        }
        else if (reciprocal && arf_is_int(mid) && arf_sgn(mid) < 0)
        {
            arb_zero(res);
            return 1;
        }
        else
        {
            /* todo: cutoffs for larger denominators */

            /* fast gamma(n), gamma(n/2) or gamma(n/4), ... */
            if (arf_cmpabs_2exp_si(mid, prec) < 0 &&
                (arf_is_int_2exp_si(mid, -2) || (prec > 1000 && arf_is_int_2exp_si(mid, -prec / 50))))
            {
                fmpq_t a;
                fmpq_init(a);
                arf_get_fmpq(a, mid);
                arb_gamma_fmpq(res, a, prec + 2 * reciprocal);
                if (reciprocal)
                    arb_inv(res, res, prec);
                fmpq_clear(a);
                return 1;
            }
        }
    }

    return 0;
}

void
arb_hypgeom_gamma_stirling(arb_t y, const arb_t x, int reciprocal, slong prec)
{
    int reflect;
    slong r, n, wp, ebits;
    arb_t t, u, v;
    double acc;

    /* for large x (if exact or accurate enough), increase precision */
    if (arf_cmpabs_2exp_si(arb_midref(x), 3) > 0)
    {
        ebits = ARF_EXP(arb_midref(x));

        if (COEFF_IS_MPZ(ebits) || ebits > 10 * prec + 4096)
        {
            arb_indeterminate(y);
            return;
        }
    }
    else
        ebits = 0;

    acc = arb_rel_accuracy_bits(x);
    acc = FLINT_MAX(acc, 0);
    wp = FLINT_MIN(prec + ebits, acc + 20);
    wp = FLINT_MAX(wp, 2);
    wp = wp + FLINT_BIT_COUNT(wp);

    if (acc < 3)  /* try to avoid divisions blowing up */
    {
        if (arf_cmp_d(arb_midref(x), -0.5) < 0)
        {
            reflect = 1;
            r = 0;
        }
        else if (arf_cmp_si(arb_midref(x), 1) < 0)
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
        arb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);
    }

    arb_init(t);
    arb_init(u);
    arb_init(v);

    if (reflect)
    {
        arb_sub_ui(t, x, 1, wp);
        arb_neg(t, t);
        arb_hypgeom_rising_ui_rec(u, t, r, wp);
        arb_const_pi(v, wp);
        arb_mul(u, u, v, wp);
        arb_add_ui(t, t, r, wp);
        arb_hypgeom_gamma_stirling_inner(v, t, n, wp);

        if (reciprocal)
        {
            /* rgamma(x) = gamma(1-x+r) sin(pi x) / ((rf(1-x, r) * pi) */
            arb_exp(v, v, wp);
            arb_sin_pi(t, x, wp);
            arb_mul(v, v, t, wp);
            arb_mul(y, u, v, wp);
            arb_div(y, v, u, prec);
        }
        else
        {
            /* gamma(x) = (rf(1-x, r) * pi) rgamma(1-x+r) csc(pi x) */
            arb_neg(v, v);
            arb_exp(v, v, wp);
            arb_csc_pi(t, x, wp);
            arb_mul(v, v, t, wp);
            arb_mul(y, v, u, prec);
        }
    }
    else
    {
        arb_add_ui(t, x, r, wp);
        arb_hypgeom_gamma_stirling_inner(u, t, n, wp);

        if (reciprocal)
        {
            /* rgamma(x) = rf(x,r) rgamma(x+r) */
            arb_neg(u, u);
            arb_exp(u, u, wp);
            arb_hypgeom_rising_ui_rec(v, x, r, wp);
            arb_mul(y, v, u, prec);
        }
        else
        {
            /* gamma(x) = gamma(x+r) / rf(x,r) */
            arb_exp(u, u, wp);
            arb_hypgeom_rising_ui_rec(v, x, r, wp);
            arb_div(y, u, v, prec);
        }
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
}

void
arb_hypgeom_gamma(arb_t y, const arb_t x, slong prec)
{
    if (arb_hypgeom_gamma_exact(y, x, 0, prec))
        return;

    if (arb_hypgeom_gamma_taylor(y, x, 0, prec))
        return;

    arb_hypgeom_gamma_stirling(y, x, 0, prec);
}

void
arb_hypgeom_rgamma(arb_t y, const arb_t x, slong prec)
{
    if (arb_hypgeom_gamma_exact(y, x, 1, prec))
        return;

    if (arb_hypgeom_gamma_taylor(y, x, 1, prec))
        return;

    arb_hypgeom_gamma_stirling(y, x, 1, prec);
}

