/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_ui_pow_ui, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t r1, r2;
        ulong a, exp;
        slong prec;

        a = n_randtest(state);
        exp = n_randtest(state);

        if (n_randint(state, 10) == 0)
            prec = 2 + n_randint(state, 5000);
        else
            prec = 2 + n_randint(state, 500);

        arb_init(r1);
        arb_init(r2);

        arb_randtest(r1, state, 1 + n_randint(state, 1000), 200);
        arb_randtest(r2, state, 1 + n_randint(state, 1000), 200);

        arb_set_ui(r1, a);
        arb_pow_ui(r1, r1, exp, prec);

        arb_ui_pow_ui(r2, a, exp, prec);

        if (!arb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = %wu\n\n", a);
            flint_printf("exp = %wu\n\n", exp);
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("r1 = "); arb_print(r1); flint_printf("\n\n");
            flint_printf("r2 = "); arb_print(r2); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(r1);
        arb_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
