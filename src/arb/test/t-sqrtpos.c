/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sqrtpos, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        slong prec = 2 + n_randint(state, 200);

        arb_init(a);
        arb_init(b);
        arb_init(c);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        } while (arb_contains_negative(a));

        arb_sqrtpos(b, a, prec);
        arb_mul(c, b, b, prec);

        if (!arb_contains(c, a) || arb_contains_negative(b))
        {
            flint_printf("FAIL: containment / positivity\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_sqrtpos(a, a, prec);

        if (!arb_equal(a, b))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
