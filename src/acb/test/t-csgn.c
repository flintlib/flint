/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_csgn, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, y;
        arb_t a;
        slong prec;

        acb_init(x);
        acb_init(y);
        arb_init(a);

        acb_randtest_special(x, state, 1 + n_randint(state, 200), 2 + n_randint(state, 100));
        arb_randtest_special(a, state, 1 + n_randint(state, 200), 2 + n_randint(state, 100));

        prec = 2 + n_randint(state, 200);

        acb_csgn(a, x);

        if (acb_is_zero(x))
        {
            acb_zero(y);
        }
        else
        {
            acb_mul(y, x, x, prec);
            acb_sqrt(y, y, prec);
            acb_div(y, y, x, prec);
        }

        if (!arb_contains(acb_realref(y), a))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("y = "); acb_printd(y, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y);
        arb_clear(a);
    }

    TEST_FUNCTION_END(state);
}
