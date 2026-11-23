/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

/* a random number with real and imaginary parts in [-scale, scale] */
static void
randc_2f1(acb_t x, flint_rand_t state, double scale, int real)
{
    acb_set_d_d(x, scale * (n_randint(state, 2001) / 1000.0 - 1.0),
        real ? 0.0 : scale * (n_randint(state, 2001) / 1000.0 - 1.0));
}

TEST_FUNCTION_START(acb_hypgeom_2f1_integration, state)
{
    slong iter, nice = 0, computed = 0;

    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, r1, r2, t;
        slong prec1, prec2;
        double scale;
        int flags, good, real;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r1); acb_init(r2); acb_init(t);

        prec1 = 20 + n_randint(state, (n_randint(state, 4) == 0) ? 300 : 100);
        prec2 = prec1 + 20 + n_randint(state, 100);
        scale = exp(log(200.0) * n_randint(state, 1001) / 1000.0);
        flags = n_randint(state, 2);
        real = (n_randint(state, 4) == 0);
        good = (n_randint(state, 5) != 0);

        /* good: Re(b) >= 1, Re(c - b) >= 1, z off [1, inf) */
        randc_2f1(a, state, 2.0 * scale, real);
        randc_2f1(b, state, scale, real);
        randc_2f1(c, state, scale, real);
        randc_2f1(z, state, 2.0, n_randint(state, 6) == 0);
        if (n_randint(state, 3) == 0)
            acb_mul_2exp_si(z, z, -2);
        if (good)
        {
            arb_abs(acb_realref(b), acb_realref(b));
            arb_abs(acb_realref(c), acb_realref(c));
            arb_add_ui(acb_realref(b), acb_realref(b), 1, 60);
            arb_add_ui(acb_realref(c), acb_realref(c), 1, 60);
            acb_add(c, c, b, 60);
            if (arb_is_zero(acb_imagref(z)) && arf_cmp_2exp_si(arb_midref(acb_realref(z)), 0) >= 0)
                acb_neg(z, z);
        }
        if (n_randint(state, 5) == 0)
        {
            arb_floor(acb_realref(a), acb_realref(a), 60);
            arb_zero(acb_imagref(a));
        }

        if (n_randint(state, 4) == 0)
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -20 - (slong) n_randint(state, 60));

        acb_hypgeom_2f1_integration(r1, a, b, c, z, flags, prec1);

        if (good)
        {
            nice++;
            computed += acb_is_finite(r1);
        }

        if (acb_is_finite(r1))
        {
            if (n_randint(state, 2))
            {
                acb_hypgeom_2f1(r2, a, b, c, z, flags, prec2 + 3 * (slong) scale + 100);
            }
            else
            {
                acb_hypgeom_2f1_integration(r2, a, b, c, z, flags, prec2);
            }

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("flags = %d, prec1 = %wd, prec2 = %wd\n\n", flags, prec1, prec2);
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r1); acb_clear(r2); acb_clear(t);
    }

    if (nice >= 20 && computed < 0.8 * nice)
    {
        flint_printf("FAIL: %wd of %wd computed\n\n", computed, nice);
        flint_abort();
    }

    /* Re(z) < 0, Re(b) >= 1 and Re(a - c) >= 0 with Re(c - b) < 1,
       where the Euler integral does not exist at 1: the two integrals
       past the point 1 / (1 - z), against the series at a precision
       which makes up for their cancellation. */
    nice = computed = 0;

    for (iter = 0; iter < 120 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, r1, r2, t;
        acb_t am, bm, cm, zm;
        slong prec1, wp, acc;
        double scale;
        int flags, real, inexact;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r1); acb_init(r2); acb_init(t);
        acb_init(am); acb_init(bm); acb_init(cm); acb_init(zm);

        prec1 = 30 + n_randint(state, (n_randint(state, 4) == 0) ? 200 : 80);
        scale = exp(log(60.0) * n_randint(state, 1001) / 1000.0);
        flags = n_randint(state, 2);
        real = (n_randint(state, 3) == 0);
        inexact = (n_randint(state, 3) == 0);

        /* b = 1 + |.|, c = b + 1 - |.| (Re(c - b) < 1), a = c + |.| */
        randc_2f1(b, state, scale, real);
        arb_abs(acb_realref(b), acb_realref(b));
        arb_add_ui(acb_realref(b), acb_realref(b), 1, 60);
        randc_2f1(t, state, scale, real);
        arb_abs(acb_realref(t), acb_realref(t));
        acb_sub(c, b, t, 60);
        arb_add_ui(acb_realref(c), acb_realref(c), 1, 60);
        if (n_randint(state, 2))
            arb_mul_2exp_si(acb_realref(c), acb_realref(c), -1);
        randc_2f1(t, state, scale, real);
        arb_abs(acb_realref(t), acb_realref(t));
        acb_add(a, c, t, 60);

        /* (a, b, c or the differences b - c, a - c integers) */
        if (n_randint(state, 4) == 0)
        {
            arb_floor(acb_realref(a), acb_realref(a), 60);
            arb_floor(acb_realref(c), acb_realref(c), 60);
            arb_ceil(acb_realref(b), acb_realref(b), 60);
            arb_zero(acb_imagref(a));
            arb_zero(acb_imagref(b));
            arb_zero(acb_imagref(c));
            acb_sub(t, a, c, 60);
            if (arb_is_negative(acb_realref(t)))
                acb_set(a, c);
        }

        randc_2f1(z, state, 4.0, real || n_randint(state, 4) == 0);
        arb_abs(acb_realref(z), acb_realref(z));
        arb_neg(acb_realref(z), acb_realref(z));
        if (n_randint(state, 3) == 0)
            acb_mul_2exp_si(z, z, -3);
        /* (strictly in the left half-plane) */
        arb_set_d(acb_realref(t), 0.001);
        arb_sub(acb_realref(z), acb_realref(z), acb_realref(t), 60);

        acb_set(am, a); acb_set(bm, b); acb_set(cm, c); acb_set(zm, z);

        if (inexact)
        {
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -20 - (slong) n_randint(state, 40));
            if (!real)
                mag_set_ui_2exp_si(arb_radref(acb_imagref(z)), 1, -20 - (slong) n_randint(state, 40));
            if (n_randint(state, 2))
            {
                mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 1, -20 - (slong) n_randint(state, 40));
                mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 1, -20 - (slong) n_randint(state, 40));
                mag_set_ui_2exp_si(arb_radref(acb_realref(c)), 1, -20 - (slong) n_randint(state, 40));
            }
        }

        acb_hypgeom_2f1_integration(r1, a, b, c, z, flags, prec1);

        /* (expected to be computed where the saddle points are complex,
           for parameters which are not small) */
        {
            _acb_hypgeom_cd A, B, C, Z;

            A.re = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
            A.im = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
            B.re = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
            B.im = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
            C.re = arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR);
            C.im = arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR);
            Z.re = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
            Z.im = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);

            if (scale >= 10.0 && !inexact && Z.re < -0.05
                && _acb_hypgeom_est_2f1_left_oscillatory(A, B, C, Z))
            {
                nice++;
                computed += acb_is_finite(r1);
            }
        }

        if (acb_is_finite(r1))
        {
            /* the series at the midpoints */
            acc = -1;
            for (wp = prec1 + 64 + 8 * (slong) scale; wp <= 10000 && acc < prec1 + 10; wp *= 2)
            {
                _acb_hypgeom_2f1_nointegration(r2, am, bm, cm, zm, flags, wp);
                acc = acb_rel_accuracy_bits(r2);
            }

            if (acc >= prec1 + 10 && !acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap (two integrals)\n\n");
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("flags = %d, prec1 = %wd\n\n", flags, prec1);
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r1); acb_clear(r2); acb_clear(t);
        acb_clear(am); acb_clear(bm); acb_clear(cm); acb_clear(zm);
    }

    if (nice >= 20 && computed < 0.6 * nice)
    {
        flint_printf("FAIL: %wd of %wd computed (two integrals)\n\n", computed, nice);
        flint_abort();
    }

    /* Re(a) <= 0 and Re(b) >= 1 with Re(c - b) < 1 inside of the circle
       |z - 1/2| = 1/2 (the functions of the type of the Jacobi
       polynomials): the loop around 1, against the series at a
       precision which makes up for their cancellation. */
    nice = computed = 0;

    for (iter = 0; iter < 120 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, r1, r2, t;
        acb_t am, bm, cm, zm;
        slong prec1, wp, acc;
        double scale, rho, phi;
        int flags, real, inexact;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r1); acb_init(r2); acb_init(t);
        acb_init(am); acb_init(bm); acb_init(cm); acb_init(zm);

        prec1 = 30 + n_randint(state, (n_randint(state, 4) == 0) ? 200 : 80);
        scale = exp(log(60.0) * n_randint(state, 1001) / 1000.0);
        flags = n_randint(state, 2);
        real = (n_randint(state, 3) == 0);
        inexact = (n_randint(state, 3) == 0);

        /* a = -|.|, b = 1 + |.|, c = b + 1 - |.| (Re(c - b) < 1) */
        randc_2f1(a, state, scale, real);
        arb_abs(acb_realref(a), acb_realref(a));
        arb_neg(acb_realref(a), acb_realref(a));
        randc_2f1(b, state, scale, real);
        arb_abs(acb_realref(b), acb_realref(b));
        arb_add_ui(acb_realref(b), acb_realref(b), 1, 60);
        randc_2f1(t, state, scale, real);
        arb_abs(acb_realref(t), acb_realref(t));
        acb_sub(c, b, t, 60);
        arb_add_ui(acb_realref(c), acb_realref(c), 1, 60);
        if (n_randint(state, 2))
            arb_mul_2exp_si(acb_realref(c), acb_realref(c), -1);

        /* (the Jacobi polynomials: a = -n, b = n + alpha + beta + 1,
           c = alpha + 1; other integers) */
        if (n_randint(state, 4) == 0)
        {
            arb_floor(acb_realref(a), acb_realref(a), 60);
            arb_zero(acb_imagref(a));
            if (n_randint(state, 2))
            {
                acb_set_d(c, 1.0 + (n_randint(state, 9) / 4.0));
                acb_sub(b, c, a, 60);
                acb_set_d(t, n_randint(state, 9) / 4.0);
                acb_add(b, b, t, 60);
            }
        }
        else if (n_randint(state, 6) == 0)
        {
            arb_floor(acb_realref(c), acb_realref(c), 60);
            arb_ceil(acb_realref(b), acb_realref(b), 60);
            arb_zero(acb_imagref(b));
            arb_zero(acb_imagref(c));
        }

        /* z = 1/2 + rho e^(i phi) / 2, rho < 1 */
        rho = 0.98 * n_randint(state, 1001) / 1000.0;
        phi = (real || n_randint(state, 4) == 0) ? 3.14159265358979323846 * n_randint(state, 2)
            : 6.28318530717958647693 * n_randint(state, 1001) / 1000.0;
        acb_set_d_d(z, 0.5 + 0.5 * rho * cos(phi), (real || sin(phi) == 0.0
            || fabs(sin(phi)) < 1e-10) ? 0.0 : 0.5 * rho * sin(phi));

        acb_set(am, a); acb_set(bm, b); acb_set(cm, c); acb_set(zm, z);

        if (inexact)
        {
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -20 - (slong) n_randint(state, 40));
            if (!real)
                mag_set_ui_2exp_si(arb_radref(acb_imagref(z)), 1, -20 - (slong) n_randint(state, 40));
            if (n_randint(state, 2))
            {
                mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 1, -20 - (slong) n_randint(state, 40));
                mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 1, -20 - (slong) n_randint(state, 40));
                mag_set_ui_2exp_si(arb_radref(acb_realref(c)), 1, -20 - (slong) n_randint(state, 40));
            }
        }

        acb_hypgeom_2f1_integration(r1, a, b, c, z, flags, prec1);

        /* (expected to be computed where the saddle points are complex,
           for parameters which are not small: the test for the two
           integrals with b - c + 1, b, a + b - c + 1 at 1 - 1/z) */
        {
            _acb_hypgeom_cd A, B, C, Z;
            double d;

            B.re = arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR);
            B.im = arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR);
            A.re = B.re - arf_get_d(arb_midref(acb_realref(c)), ARF_RND_NEAR) + 1.0;
            A.im = B.im - arf_get_d(arb_midref(acb_imagref(c)), ARF_RND_NEAR);
            C.re = A.re + arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
            C.im = A.im + arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
            Z.re = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
            Z.im = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
            d = Z.re * Z.re + Z.im * Z.im;
            Z.re = 1.0 - Z.re / d;
            Z.im = Z.im / d;

            if (scale >= 10.0 && !inexact && Z.re < -0.05
                && _acb_hypgeom_est_2f1_left_oscillatory(A, B, C, Z))
            {
                nice++;
                computed += acb_is_finite(r1);
            }
        }

        if (acb_is_finite(r1))
        {
            /* the series at the midpoints */
            acc = -1;
            for (wp = prec1 + 64 + 8 * (slong) scale; wp <= 10000 && acc < prec1 + 10; wp *= 2)
            {
                _acb_hypgeom_2f1_nointegration(r2, am, bm, cm, zm, flags, wp);
                acc = acb_rel_accuracy_bits(r2);
            }

            if (acc >= prec1 + 10 && !acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap (loop)\n\n");
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("flags = %d, prec1 = %wd\n\n", flags, prec1);
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r1); acb_clear(r2); acb_clear(t);
        acb_clear(am); acb_clear(bm); acb_clear(cm); acb_clear(zm);
    }

    if (nice >= 20 && computed < 0.6 * nice)
    {
        flint_printf("FAIL: %wd of %wd computed (loop)\n\n", computed, nice);
        flint_abort();
    }

    /* Re(b) >= 1 and Re(c - b) >= 1 for z on the cut (1, inf), where
       the function is the limit from the lower half-plane, and close
       to it on both sides: the Euler integral over a polygon past the
       point 1 / z, against the series at the midpoints. */
    nice = computed = 0;

    for (iter = 0; iter < 120 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, r1, r2, t;
        acb_t am, bm, cm, zm;
        slong prec1, wp, acc;
        double scale;
        int flags, real, inexact, where;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r1); acb_init(r2); acb_init(t);
        acb_init(am); acb_init(bm); acb_init(cm); acb_init(zm);

        prec1 = 30 + n_randint(state, (n_randint(state, 4) == 0) ? 200 : 80);
        scale = exp(log(60.0) * n_randint(state, 1001) / 1000.0);
        flags = n_randint(state, 2);
        real = (n_randint(state, 3) == 0);
        inexact = (n_randint(state, 3) == 0);
        where = n_randint(state, 4);

        /* b = 1 + |.|, c = b + 1 + |.| */
        randc_2f1(a, state, 2.0 * scale, real);
        randc_2f1(b, state, scale, real);
        arb_abs(acb_realref(b), acb_realref(b));
        arb_add_ui(acb_realref(b), acb_realref(b), 1, 60);
        randc_2f1(c, state, scale, real);
        arb_abs(acb_realref(c), acb_realref(c));
        arb_add_ui(acb_realref(c), acb_realref(c), 1, 60);
        acb_add(c, c, b, 60);

        if (n_randint(state, 4) == 0)
        {
            arb_floor(acb_realref(a), acb_realref(a), 60);
            arb_floor(acb_realref(b), acb_realref(b), 60);
            arb_floor(acb_realref(c), acb_realref(c), 60);
            arb_zero(acb_imagref(a));
            arb_zero(acb_imagref(b));
            arb_zero(acb_imagref(c));
        }

        /* z = 1 + |.| on the cut, or off it by 2^-200 or by 2^-10 */
        acb_set_d(z, 1.001 + 3.0 * n_randint(state, 1001) / 1000.0);
        if (n_randint(state, 3) == 0)
            acb_set_d(z, 1.001 + 0.1 * n_randint(state, 1001) / 1000.0);
        if (where != 0)
        {
            arb_one(acb_imagref(z));
            arb_mul_2exp_si(acb_imagref(z), acb_imagref(z), (where == 3) ? -10 : -200);
            if (n_randint(state, 2))
                arb_neg(acb_imagref(z), acb_imagref(z));
        }

        acb_set(am, a); acb_set(bm, b); acb_set(cm, c); acb_set(zm, z);

        if (inexact)
        {
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -20 - (slong) n_randint(state, 40));
            if (n_randint(state, 2))
            {
                mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 1, -20 - (slong) n_randint(state, 40));
                mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 1, -20 - (slong) n_randint(state, 40));
                mag_set_ui_2exp_si(arb_radref(acb_realref(c)), 1, -20 - (slong) n_randint(state, 40));
            }
        }

        acb_hypgeom_2f1_integration(r1, a, b, c, z, flags, prec1);

        /* (expected to be computed for parameters which are not small:
           the polygon itself is not integrated) */
        if (scale >= 10.0 && !inexact)
        {
            nice++;
            computed += acb_is_finite(r1);
        }

        if (acb_is_finite(r1))
        {
            /* the series at the midpoints */
            acc = -1;
            for (wp = prec1 + 64 + 8 * (slong) scale; wp <= 10000 && acc < prec1 + 10; wp *= 2)
            {
                _acb_hypgeom_2f1_nointegration(r2, am, bm, cm, zm, flags, wp);
                acc = acb_rel_accuracy_bits(r2);
            }

            if (acc >= prec1 + 10 && !acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap (cut)\n\n");
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("flags = %d, prec1 = %wd\n\n", flags, prec1);
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r1); acb_clear(r2); acb_clear(t);
        acb_clear(am); acb_clear(bm); acb_clear(cm); acb_clear(zm);
    }

    if (nice >= 20 && computed < 0.6 * nice)
    {
        flint_printf("FAIL: %wd of %wd computed (cut)\n\n", computed, nice);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
