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

/* z, tau, zeta(z, tau) checked with Mathematica:
   N[{z, tau, WeierstrassZeta[z, WeierstrassInvariants[{1, tau}/2]]}, 20] */

/* testdata is defined in t-sigma.c and t-zeta.c */
#define testdata testdata1
const double testdata[NUM_TESTS][6] = {
    { 1.4142135623730950488, 1.7320508075688772935,
      2.2360679774997896964, 2.6457513110645905905,
      4.6708530542187465956, 2.5654839019069906619 },
    { -3.0, -2.0, -7.0, 3.0, -9.8696028584909729616, -3.4498761151331612449 },
    { 0.1, 0.0, 0.6, 0.2, 10.0222822996630537305, -0.0754856899565257555 },
    { 0.0, 0.1, 0.6, 0.2, -0.0754856899565257555, -10.0222822996630537305 },
    { 0.5, 0.0, 0.333333333333333333, 20.0, 1.6449340668482264365, 0.0 },
    { 3.6666666666666667, 2014.0, -3.1415926535897932385, 0.1,
      -147965.49261828128352, -10426.83705290890375 }
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

TEST_FUNCTION_START(acb_elliptic_zeta, state)
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
            if (i == NUM_TESTS - 1) /* sensitive to rounding errors in doubles */
                acb_set_dddd(p2, testdata[i][4], 1e-6, testdata[i][5], 1e-6);
            else
                acb_set_dddd(p2, testdata[i][4], EPS, testdata[i][5], EPS);

            acb_elliptic_zeta(p1, z, tau, 2 + n_randint(state, 400));

            if (!acb_overlaps(p1, p2))
            {
                flint_printf("FAIL (test value)\n");
                flint_printf("tau = "); acb_printd(tau, 15); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 15); flint_printf("\n\n");
                flint_printf("p1 = "); acb_printd(p1, 15); flint_printf("\n\n");
                flint_printf("p2 = "); acb_printd(p2, 15); flint_printf("\n\n");
                flint_abort();
            }

            acb_elliptic_zeta(p2, z, tau, 2 + n_randint(state, 800));

            if (!acb_overlaps(p1, p2))
            {
                flint_printf("FAIL (test value 2)\n");
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
        acb_t tau, z1, z2, e1, e2, p1, p2;
        slong m, n, e0, prec0, prec1, prec2;

        acb_init(tau);
        acb_init(z1);
        acb_init(z2);
        acb_init(e1);
        acb_init(e2);
        acb_init(p1);
        acb_init(p2);

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
        acb_mul_2exp_si(e1, e1, 1);
        acb_mul_2exp_si(e2, e2, 1);

        acb_randtest(z1, state, prec0, e0);
        acb_randtest(p1, state, prec0, e0);
        acb_randtest(p2, state, prec0, e0);

        /* z2 = z1 + m + n*tau */
        m = n_randint(state, 10);
        n = n_randint(state, 10);
        acb_add_ui(z2, z1, m, prec0);
        acb_addmul_ui(z2, tau, n, prec0);

        /* zeta(z + 1) = zeta(z) + 2 zeta(1) */
        /* zeta(z + tau) = zeta(z) + 2 zeta(tau) */

        acb_elliptic_zeta(p1, z1, tau, prec1);
        acb_elliptic_zeta(p2, z2, tau, prec2);

        acb_submul_ui(p2, e1, m, prec2);
        acb_submul_ui(p2, e2, n, prec2);

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
    }

    TEST_FUNCTION_END(state);
}
#undef testdata1
