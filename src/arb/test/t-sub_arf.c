/*
    Copyright (C) 2012-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sub_arf, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        arf_t x;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arf_init(x);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        arf_randtest_special(x, state, 1 + n_randint(state, 2000), 100);

        prec = 2 + n_randint(state, 2000);

        arb_set_arf(b, x);
        arb_sub_arf(c, a, x, prec);
        arb_sub(d, a, b, prec);

        if (!arb_equal(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arf_clear(x);
    }

    /* aliasing */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        arf_t x;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arf_init(x);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        arf_randtest_special(x, state, 1 + n_randint(state, 2000), 100);

        prec = 2 + n_randint(state, 2000);

        arb_set_arf(b, x);
        arb_sub_arf(c, a, x, prec);
        arb_sub_arf(a, a, x, prec);

        if (!arb_equal(a, c))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arf_clear(x);
    }

    TEST_FUNCTION_END(state);
}
