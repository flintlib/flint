/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

void
acb_dirichlet_lerch_phi_test(acb_t res, const acb_t z, const acb_t s, const acb_t a, int algorithm, slong prec)
{
    switch (algorithm % 3)
    {
        case 0:
            acb_dirichlet_lerch_phi_direct(res, z, s, a, prec);
            if (!acb_is_finite(res))
                acb_dirichlet_lerch_phi(res, z, s, a, prec);
            break;
        case 1:
            acb_dirichlet_lerch_phi_integral(res, z, s, a, prec);
            break;
        default:
            acb_dirichlet_lerch_phi(res, z, s, a, prec);
            break;
    }
}

static void
acb_set_dddd(acb_t z, double a, double ar, double b, double br)
{
    arf_set_d(arb_midref(acb_realref(z)), a);
    mag_set_d(arb_radref(acb_realref(z)), ar);
    arf_set_d(arb_midref(acb_imagref(z)), b);
    mag_set_d(arb_radref(acb_imagref(z)), br);
}

#define NUM_TESTS 15

/* z, s, a, phi(z, s, a) */
const double testdata[NUM_TESTS][8] = {
    { 1.0, 0.00390625, 1.5, 0.0, 2.25, 0.0, 1.3422063426254739626, 0.13465412214420029791 },
    { 1.0, 0.0, 1.5, 0.0, 2.25, 0.0, 1.4975136076666680360, 0.0 },
    { 1.0, -0.00390625, 1.5, 0.0, 2.25, 0.0, 1.3422063426254739626, -0.13465412214420029791 },
    { 2.0, 0.00390625, 1.5, 0.0, 2.25, 0.0, -0.16955701974487404975, 0.61946301461898363407 },
    { 2.0, 0.0, 1.5, 0.0, 2.25, 0.0, -0.17140382129993246416, -0.62044054732729604008 },
    { 2.0, -0.00390625, 1.5, 0.0, 2.25, 0.0, -0.16955701974487404975, -0.61946301461898363407 },
    { 0.6875, 0.0, 0.0, 0.0, -1.0, 2.0, 3.2, 0.0 },
    { 0.0, 8.0, 1.0, 0.0, 1.0, -1.0, -0.21201062033942531891, 0.15142847985530888328 },
    { 0.0, 8.0, 1.0, 1.0, 1.0, -1.0, -0.18714764709994647918, -0.031327583631588241802 },
    { 1.875, 2.625, -1.0, 0.0, 1.0, 0.0, -0.10448979591836734694, -0.078367346938775510204 },
    { -2.5, -1.0, -3.75, 0.0, -3.75, 7.5, 195.9465378716877103, 889.34175804326870925 },
    { 0.0, -1.5, 1.0, 0.5, -1.0, 6.0, -0.20452787205323008395, 0.14062505401787546806 },
    { -3.0, 3.0, 9.5, 0.0, 0.0, -5.5, 5.5775511785583065106e-7, -4.488380385476415194e-7 },
    { -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -0.75, 0.0 },
    { -3.0, 0.0, -2.0, 0.0, -2.0, 0.0, 1.84375, 0.0 },
};

TEST_FUNCTION_START(acb_dirichlet_lerch_phi, state)
{
    slong iter;

    /* check test values */
    for (iter = 0; iter < 5 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong i;
        acb_t z, s, a, p1, p2;
        int alg;

        acb_init(z);
        acb_init(s);
        acb_init(a);
        acb_init(p1);
        acb_init(p2);

        for (i = 0; i < NUM_TESTS; i++)
        {
            alg = n_randlimb(state);

            acb_set_dddd(z, testdata[i][0], 0.0, testdata[i][1], 0.0);
            acb_set_dddd(s, testdata[i][2], 0.0, testdata[i][3], 0.0);
            acb_set_dddd(a, testdata[i][4], 1e-14, testdata[i][5], 1e-14);
            acb_set_dddd(p2, testdata[i][6], 1e-14, testdata[i][7], 1e-14);

            acb_dirichlet_lerch_phi_test(p1, z, s, a, alg, 2 + n_randint(state, 100));

            if (!acb_overlaps(p1, p2))
            {
                flint_printf("FAIL (test value)\n");
                flint_printf("z = "); acb_printd(z, 15); flint_printf("\n\n");
                flint_printf("s = "); acb_printd(s, 15); flint_printf("\n\n");
                flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
                flint_printf("p1 = "); acb_printd(p1, 15); flint_printf("\n\n");
                flint_printf("p2 = "); acb_printd(p2, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(z);
        acb_clear(s);
        acb_clear(a);
        acb_clear(p1);
        acb_clear(p2);
    }

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z, s, a, b, r1, r2, r3;
        slong prec1, prec2;
        int alg1, alg2, alg3;

        /* flint_printf("iter %wd\n", iter); */

        acb_init(z);
        acb_init(s);
        acb_init(a);
        acb_init(b);
        acb_init(r1);
        acb_init(r2);
        acb_init(r3);

        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);
        alg1 = n_randlimb(state);
        alg2 = n_randlimb(state);
        alg3 = n_randlimb(state);

        acb_randtest(z, state, 1 + n_randint(state, 500), 3);
        acb_randtest(s, state, 1 + n_randint(state, 500), 3);
        acb_randtest(a, state, 1 + n_randint(state, 500), 3);
        acb_randtest(b, state, 1 + n_randint(state, 500), 10);
        acb_randtest(r1, state, 1 + n_randint(state, 500), 10);
        acb_randtest(r2, state, 1 + n_randint(state, 500), 10);

        if (n_randint(state, 2))
            acb_set_si(z, n_randint(state, 5) - 2);
        if (n_randint(state, 2))
            acb_set_si(a, n_randint(state, 5) - 2);
        if (n_randint(state, 2))
            acb_set_si(s, n_randint(state, 5) - 2);

        acb_dirichlet_lerch_phi_test(r1, z, s, a, alg1, prec1);
        acb_dirichlet_lerch_phi_test(r2, z, s, a, alg2, prec2);

        if (!acb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd\n\n", iter);
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        /* test phi(z,s,a) = z*phi(z,s,a+1) + a^-s */
        acb_add_ui(b, a, 1, prec2);
        acb_dirichlet_lerch_phi_test(r2, z, s, b, alg3, prec2);
        acb_neg(r3, s);
        acb_pow(r3, a, r3, prec2);
        acb_addmul(r3, r2, z, prec2);

        if (!acb_overlaps(r1, r3))
        {
            flint_printf("FAIL (2): overlap\n\n");
            flint_printf("iter = %wd\n\n", iter);
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_printf("r3 = "); acb_printd(r3, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z);
        acb_clear(s);
        acb_clear(a);
        acb_clear(b);
        acb_clear(r1);
        acb_clear(r2);
        acb_clear(r3);
    }

    TEST_FUNCTION_END(state);
}
