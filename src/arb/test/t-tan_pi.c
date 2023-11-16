/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_tan_pi, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(x);
        arb_init(y);
        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        arb_tan_pi(a, x, prec1);
        arb_tan_pi(b, x, prec2);

        /* check consistency */
        if (!arb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        /* compare with tan */
        arb_const_pi(c, prec1);
        arb_mul(y, x, c, prec1);
        arb_tan(c, y, prec1);

        if (!arb_overlaps(a, c))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_tan_pi(x, x, prec1);

        if (!arb_overlaps(a, x))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    /* Check a few special points. */
    {
        arb_t x, y;
        slong i, prec;

        arb_init(x);
        arb_init(y);

        prec = 32;

        for (i = -1; i <= 1; i++)
        {
            int arbitrary_integer = 7;

            arb_set_d(x, arbitrary_integer*i);
            arb_tan_pi(y, x, prec);
            if (!arb_is_zero(y))
            {
                flint_printf("FAIL: (integer)\n");
                flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
                flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
                flint_abort();
            }

            arb_set_d(x, arbitrary_integer*i + 0.5);
            arb_tan_pi(y, x, prec);
            if (arb_is_finite(y))
            {
                flint_printf("FAIL: (integer + 1/2)\n");
                flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
                flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
                flint_abort();
            }
        }
    }

    TEST_FUNCTION_END(state);
}
