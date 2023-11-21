/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sqrt1pm1, state)
{
    slong iter;

    for (iter = 0; iter < 20000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        slong prec0, prec1, prec2;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 1000;

        prec1 = 2 + n_randint(state, prec0);
        prec2 = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, prec0), 100);

        arb_sqrt1pm1(b, a, prec1);
        arb_sqrt1pm1(c, a, prec2);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        /* compare with sqrt */
        arb_add_ui(d, a, 1, prec2);
        arb_sqrt(d, d, prec2);
        arb_sub_ui(d, d, 1, prec2);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: comparison with log\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_sqrt1pm1(a, a, prec1);

        if (!arb_overlaps(a, b))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
