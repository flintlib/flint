/*
    Copyright (C) 2015 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sinc, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, a;
        slong prec, mag_bits;

        prec = 2 + n_randint(state, 500);
        mag_bits = 6 + n_randint(state, 500);

        arb_init(x);
        arb_randtest(x, state, 1 + n_randint(state, prec), mag_bits);

        arb_init(a);
        arb_sinc(a, x, prec);

        /* sinc(x) approximately in [-2, 2] */
        {
            arb_t b;
            arb_init(b);
            arb_zero(b);
            mag_set_d(arb_radref(b), 1.01);
            mag_mul_2exp_si(arb_radref(b), arb_radref(b), 1);

            if (!arb_contains(b, a))
            {
                flint_printf("FAIL: containment\n\n");
                flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                flint_printf("a = "); arb_print(a); flint_printf("\n\n");
                flint_printf("b = "); arb_print(b); flint_printf("\n\n");
                flint_abort();
            }

            arb_clear(b);
        }

        /* sinc(x) * x = sin(x) */
        {
            arb_t b, c;
            arb_init(b);
            arb_init(c);

            arb_mul(b, a, x, prec);
            arb_sin(c, x, prec);

            if (!arb_overlaps(b, c))
            {
                flint_printf("FAIL: overlap (sin)\n\n");
                flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                flint_printf("b = "); arb_print(b); flint_printf("\n\n");
                flint_printf("c = "); arb_print(c); flint_printf("\n\n");
                flint_abort();
            }

            arb_clear(b);
            arb_clear(c);
        }

        /* aliasing */
        {
            arb_t b;
            arb_init(b);
            arb_set(b, x);
            arb_sinc(b, b, prec);
            if (!arb_equal(a, b))
            {
                flint_printf("FAIL: aliasing\n\n");
                flint_abort();
            }
            arb_clear(b);
        }

        arb_clear(x);
        arb_clear(a);
    }

    /* Check a few special intervals. */
    {
        arb_t x, y;
        slong prec;

        arb_init(x);
        arb_init(y);

        prec = 32;

        arb_neg_inf(x);
        arb_sinc(y, x, prec);
        if (!arb_is_zero(y))
        {
            flint_printf("FAIL: negative infinity\n");
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_pos_inf(x);
        arb_sinc(y, x, prec);
        if (!arb_is_zero(y))
        {
            flint_printf("FAIL: positive infinity\n");
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_zero_pm_inf(x);
        arb_sinc(y, x, prec);
        if (!(arb_is_finite(y) &&
              arb_contains_negative(y) &&
              arb_contains_positive(y) &&
              arb_contains_zero(y)))
        {
            flint_printf("FAIL: the whole extended real line\n");
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
    }

    TEST_FUNCTION_END(state);
}
