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

TEST_FUNCTION_START(arb_atanh, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, a, b;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(x);
        arb_init(a);
        arb_init(b);

        arb_randtest_special(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        arb_randtest_special(a, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        arb_atanh(a, x, prec1);
        arb_atanh(b, x, prec2);

        /* check consistency */
        if (!arb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); arb_printd(x, 15); flint_printf("\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
            flint_abort();
        }

        /* check tanh(atanh(x)) = x */
        arb_tanh(b, b, prec1);

        if (!arb_contains(b, x))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("x = "); arb_printd(x, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_atanh(x, x, prec1);

        if (!arb_overlaps(a, x))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("x = "); arb_printd(x, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(a);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
