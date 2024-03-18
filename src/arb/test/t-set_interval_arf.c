/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_set_interval_arf, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x;
        arf_t a, b;

        arb_init(x);
        arf_init(a);
        arf_init(b);

        arf_randtest_special(a, state, 200, 10);
        arf_randtest_special(b, state, 200, 10);
        if (arf_cmp(a, b) > 0)
            arf_swap(a, b);

        arb_set_interval_arf(x, a, b, 2 + n_randint(state, 200));

        if ((!arb_contains_arf(x, a) || !arb_contains_arf(x, b)) ||
            (!arf_is_nan(a) && !arf_is_nan(b) && arf_is_nan(arb_midref(x))))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arf_clear(a);
        arf_clear(b);
    }

    TEST_FUNCTION_END(state);
}
