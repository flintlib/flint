/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

/* a random exact number of the size 10^e: an integer or with a few
   fractional bits */
static void
_randbig_u(arb_t x, flint_rand_t state, slong e)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(x, 1 + n_randint(state, 1000000));
    arb_mul_2exp_si(x, x, -10);
    arb_set_ui(t, 10);
    arb_pow_ui(t, t, e, 4000);
    arb_mul(x, x, t, 4000);
    if (n_randint(state, 2))
        arb_floor(x, x, 4000);
    arb_clear(t);
}

TEST_FUNCTION_START(arb_hypgeom_u_integration, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, z, r1, r2;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 80);
        prec2 = 2 + n_randint(state, 150);

        arb_init(a);
        arb_init(b);
        arb_init(z);
        arb_init(r1);
        arb_init(r2);

        arb_randtest_precise(a, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest_precise(b, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest_precise(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest(r1, state, 1 + n_randint(state, 200), 1 + n_randint(state, 5));
        arb_randtest(r2, state, 1 + n_randint(state, 200), 1 + n_randint(state, 5));

        if (n_randint(state, 2))
            arb_add_ui(a, a, n_randint(state, 100), prec1 + 100);

        arb_add_ui(b, b, n_randint(state, 200), prec1 + 100);
        arb_add_ui(z, z, n_randint(state, 100), prec1 + 100);

        arb_hypgeom_u_integration(r1, a, b, z, prec1);

        if (arb_is_finite(r1))
        {
            if (n_randint(state, 2))
                arb_hypgeom_u(r2, a, b, z, prec2);
            else
                arb_hypgeom_u_integration(r2, a, b, z, prec2);

            if (!arb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 30); flint_printf("\n\n");
                flint_printf("z = "); arb_printd(z, 30); flint_printf("\n\n");
                flint_printf("r1 = "); arb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (iter == 0)
        {
            prec1 = 333;

            arb_set_str(a, "1.2e7", prec1);
            arb_set_str(b, "1.3e8", prec1);
            arb_set_str(z, "1.4e8", prec1);

            arb_hypgeom_u_integration(r1, a, b, z, prec1);
            arb_set_str(r2, "7.423720668435352497136323553980928479792327057967344627374385949694801879462740896706198120e-90723524 +/- 6.47e-90723615", prec1);

            if (!arb_overlaps(r1, r2) || arb_rel_accuracy_bits(r1) < arb_rel_accuracy_bits(r2) - 10)
            {
                flint_printf("FAIL: overlap (1)\n\n");
                flint_printf("r1 = "); arb_printd(r1, 100); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printd(r2, 100); flint_printf("\n\n");
                flint_abort();
            }

            arb_set_str(a, "1.0", prec1);
            arb_set_str(b, "1351.25", prec1);
            arb_set_str(z, "5118.125", prec1);

            arb_hypgeom_u_integration(r1, a, b, z, prec1);
            arb_set_str(r2, "0.0002653059837025369361090415756761601958818507141969057844945316037867908521084665825645386229229474 +/- 5.74e-101", prec1);

            if (!arb_overlaps(r1, r2) || arb_rel_accuracy_bits(r1) < arb_rel_accuracy_bits(r2) - 10)
            {
                flint_printf("FAIL: overlap (2)\n\n");
                flint_printf("r1 = "); arb_printd(r1, 100); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printd(r2, 100); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(z);
        arb_clear(r1);
        arb_clear(r2);
    }

    /* Large parameters: Kummer's transformation
       U(a, b, z) = z^(1 - b) U(a - b + 1, 2 - b, z), and U(a, a + 1, z) = z^(-a) */
    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, z, r1, r2, t, u;
        slong prec1, prec2, e;
        int which;

        arb_init(a); arb_init(b); arb_init(z); arb_init(r1); arb_init(r2); arb_init(t); arb_init(u);

        prec1 = 30 + n_randint(state, 200);
        prec2 = prec1 + 10 + n_randint(state, 100);
        e = 3 + n_randint(state, n_randint(state, 4) ? 40 : 110);
        which = n_randint(state, 4) == 0;

        _randbig_u(a, state, e - n_randint(state, 3));
        _randbig_u(b, state, e - n_randint(state, 3));
        _randbig_u(z, state, e - n_randint(state, 3));
        arb_add_ui(a, a, 1, 4000);

        if (which == 0)
        {
            /* (both sides by the integral: b <= a) */
            arb_sub(b, a, b, 4000);
            arb_hypgeom_u_integration(r1, a, b, z, prec1);
            arb_sub(t, a, b, 4000);
            arb_add_ui(t, t, 1, 4000);
            arb_sub_ui(u, b, 2, 4000);
            arb_neg(u, u);
            arb_hypgeom_u_integration(r2, t, u, z, prec2);
            arb_sub_ui(u, b, 1, 4000);
            arb_neg(u, u);
            arb_pow(t, z, u, prec2 + 4000);
            arb_mul(r2, r2, t, prec2 + 400);
        }
        else
        {
            arb_add_ui(b, a, 1, 4000);
            arb_hypgeom_u_integration(r1, a, b, z, prec1);
            arb_neg(t, a);
            arb_pow(r2, z, t, prec2 + 4000);
        }

            if (!arb_is_finite(r1) || !arb_is_finite(r2) || !arb_overlaps(r1, r2)
                || arb_rel_accuracy_bits(r1) < prec1 - 15)
            {
                flint_printf("FAIL: large parameters (%d)\n\n", which);
                flint_printf("a = "); arb_printd(a, 50); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 50); flint_printf("\n\n");
                flint_printf("z = "); arb_printd(z, 50); flint_printf("\n\n");
                flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
                flint_printf("r1 = "); arb_printd(r1, 50); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printd(r2, 50); flint_printf("\n\n");
                flint_abort();
            }

        arb_clear(a); arb_clear(b); arb_clear(z); arb_clear(r1); arb_clear(r2); arb_clear(t); arb_clear(u);
    }

    /* a high precision, where the integrand in the outer cells is below
       the range of doubles relative to its maximum */
    {
        arb_t a, b, z, r1, r2;
        arb_init(a); arb_init(b); arb_init(z); arb_init(r1); arb_init(r2);
        arb_set_d(a, 1000.25);
        arb_set_d(b, 3000.75);
        arb_set_d(z, 2000.125);
        arb_hypgeom_u_integration(r1, a, b, z, 2500);
        arb_hypgeom_u_integration(r2, a, b, z, 1200);
        if (!arb_overlaps(r1, r2) || arb_rel_accuracy_bits(r1) < 2470 || arb_rel_accuracy_bits(r2) < 1170)
        {
            flint_printf("FAIL: high precision\n\n");
            flint_printf("accuracy: %wd, %wd\n\n", arb_rel_accuracy_bits(r1), arb_rel_accuracy_bits(r2));
            flint_abort();
        }
        arb_clear(a); arb_clear(b); arb_clear(z); arb_clear(r1); arb_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
