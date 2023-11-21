/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sgn, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        int result;

        arb_init(a);
        arb_init(b);

        arb_randtest_special(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest_special(b, state, 1 + n_randint(state, 200), 10);
        arb_sgn(b, a);

        result = 1;
        if (arb_contains_zero(a))
            result = result & arb_contains_si(b, 0);
        if (arb_contains_positive(a))
            result = result & arb_contains_si(b, 1);
        if (arb_contains_negative(a))
            result = result & arb_contains_si(b, -1);

        if (!result)
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
