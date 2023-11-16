/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_modular.h"
#include "acb_elliptic.h"

#define NUM_TESTS 6
#define EPS 1e-13

/* z, tau, sigma(z, tau) checked with Mathematica:
   N[{z, tau, WeierstrassSigma[z, WeierstrassInvariants[{1, tau}/2]]}, 20] */

/* testdata is defined in t-sigma.c and t-zeta.c */
#define testdata testdata0
const double testdata[NUM_TESTS][6] = {
    { 1.4142135623730950488, 1.7320508075688772935,
      2.2360679774997896964, 2.6457513110645905905,
      3.2497809387982239642, -6.2896427497987532326 },
    { -1.0, -2.0, -1.0, 3.0, -0.17877885105742438172, 0.58508579024326876042 },
    { 0.1, 0.0, 0.6, 0.2, 0.100055263033144515447, -0.000188998253739903104 },
    { 0.0, 0.2, 0.2, 0.1, -0.05083547794781899013, 0.19409530512485630787 },
    { 0.5, 0.0, 0.333333333333333333, 20.0, 0.48022700051193809297, 0.0 },
    { 0.6666666666666667, 1.0, -3.1415926535897932385, 1.0,
      1.3092907491438550394, 0.9063920053572463817 }
};

/* acb_set_dddd is defined in t-sigma.c and t-zeta.c */
#ifndef acb_set_dddd
#define acb_set_dddd acb_set_dddd
static void
acb_set_dddd(acb_t z, double a, double ar, double b, double br)
{
    arf_set_d(arb_midref(acb_realref(z)), a);
    mag_set_d(arb_radref(acb_realref(z)), ar);
    arf_set_d(arb_midref(acb_imagref(z)), b);
    mag_set_d(arb_radref(acb_imagref(z)), br);
}
#endif

TEST_FUNCTION_START(acb_elliptic_sigma, state)
{
    slong iter;

    /* check test values */
    for (iter = 0; iter < 20 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong i;

        acb_t z, tau, p1, p2;

        acb_init(z);
        acb_init(tau);
        acb_init(p1);
        acb_init(p2);

        for (i = 0; i < NUM_TESTS; i++)
        {
            acb_set_dddd(z, testdata[i][0], 0.0, testdata[i][1], 0.0);
            acb_set_dddd(tau, testdata[i][2], 0.0, testdata[i][3], 0.0);
            acb_set_dddd(p2, testdata[i][4], EPS, testdata[i][5], EPS);

            acb_elliptic_sigma(p1, z, tau, 2 + n_randint(state, 400));

            if (!acb_overlaps(p1, p2))
            {
                flint_printf("FAIL (test value)\n");
                flint_printf("tau = "); acb_printd(tau, 15); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 15); flint_printf("\n\n");
                flint_printf("p1 = "); acb_printd(p1, 15); flint_printf("\n\n");
                flint_printf("p2 = "); acb_printd(p2, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(z);
        acb_clear(tau);
        acb_clear(p1);
        acb_clear(p2);
    }

    /* Test periods */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t tau, z1, z2, e1, e2, p1, p2, t, u;
        slong m, n, e0, prec0, prec1, prec2;

        acb_init(tau);
        acb_init(z1);
        acb_init(z2);
        acb_init(e1);
        acb_init(e2);
        acb_init(p1);
        acb_init(p2);
        acb_init(t);
        acb_init(u);

        e0 = 1 + n_randint(state, 10);
        prec0 = 2 + n_randint(state, 400);
        prec1 = 2 + n_randint(state, 400);
        prec2 = 2 + n_randint(state, 400);

        acb_randtest(tau, state, prec0, e0);
        if (arf_sgn(arb_midref(acb_imagref(tau))) < 0)
            acb_neg(tau, tau);

        acb_one(e1);
        acb_mul_2exp_si(e1, e1, -1);
        acb_elliptic_zeta(e1, e1, tau, prec0);
        acb_mul_2exp_si(e2, tau, -1);
        acb_elliptic_zeta(e2, e2, tau, prec0);

        acb_randtest(z1, state, prec0, e0);
        acb_randtest(p1, state, prec0, e0);
        acb_randtest(p2, state, prec0, e0);

        /* z2 = z1 + m + n*tau */
        m = n_randint(state, 10);
        n = n_randint(state, 10);
        acb_add_ui(z2, z1, m, prec0);
        acb_addmul_ui(z2, tau, n, prec0);

        /* sigma(z + m + n*tau) = (-1)^(m+n+mn) exp((m e1 + n e2)(m + n tau + 2z) sigma(z) */
        acb_elliptic_sigma(p1, z1, tau, prec1);
        acb_elliptic_sigma(p2, z2, tau, prec2);

        acb_mul_ui(t, e1, m, prec2);
        acb_addmul_ui(t, e2, n, prec2);
        acb_mul_2exp_si(u, z1, 1);
        acb_addmul_ui(u, tau, n, prec2);
        acb_add_ui(u, u, m, prec2);
        acb_mul(t, t, u, prec2);
        acb_neg(t, t);
        acb_exp(t, t, prec2);
        if ((m + n + m*n) % 2 == 1)
            acb_neg(t, t);
        acb_mul(p2, p2, t, prec2);

        if (!acb_overlaps(p1, p2))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("tau = "); acb_printd(tau, 15); flint_printf("\n\n");
            flint_printf("z1 = "); acb_printd(z1, 15); flint_printf("\n\n");
            flint_printf("z2 = "); acb_printd(z2, 15); flint_printf("\n\n");
            flint_printf("p1 = "); acb_printd(p1, 15); flint_printf("\n\n");
            flint_printf("p2 = "); acb_printd(p2, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau);
        acb_clear(z1);
        acb_clear(z2);
        acb_clear(e1);
        acb_clear(e2);
        acb_clear(p1);
        acb_clear(p2);
        acb_clear(t);
        acb_clear(u);
    }

    TEST_FUNCTION_END(state);
}
#undef testdata
