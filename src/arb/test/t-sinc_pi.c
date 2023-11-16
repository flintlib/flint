/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sinc_pi, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z;
        slong prec;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        prec = 2 + n_randint(state, 200);

        arb_randtest(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arb_randtest(y, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        if (n_randint(state, 2))
        {
            arb_set(y, x);
            arb_sinc_pi(y, y, prec);
        }
        else
            arb_sinc_pi(y, x, prec);

        arb_const_pi(z, prec);
        arb_mul(z, z, x, prec);
        arb_sinc(z, z, prec);

        if (!arb_overlaps(y, z))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 30); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_is_int(x))
        {
            if (arb_is_zero(x) && !arb_is_one(y))
            {
                flint_printf("FAIL: zero\n");
                flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
                flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
                flint_abort();
            }
            if (arb_is_positive(x) && !arb_is_zero(y))
            {
                flint_printf("FAIL: positive integer\n");
                flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
                flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
                flint_abort();
            }
            if (arb_is_negative(x) && !arb_is_zero(y))
            {
                flint_printf("FAIL: negative integer\n");
                flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
                flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    /* Check a few special intervals. */
    {
        arb_t x, y;
        slong prec;

        arb_init(x);
        arb_init(y);

        prec = 32;

        arb_neg_inf(x);
        arb_sinc_pi(y, x, prec);
        if (!arb_is_zero(y))
        {
            flint_printf("FAIL: negative infinity\n");
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_pos_inf(x);
        arb_sinc_pi(y, x, prec);
        if (!arb_is_zero(y))
        {
            flint_printf("FAIL: positive infinity\n");
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_zero_pm_inf(x);
        arb_sinc_pi(y, x, prec);
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
