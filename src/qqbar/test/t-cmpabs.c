/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_cmpabs, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, xr, yr;
        int c1, c2;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(xr);
        qqbar_init(yr);

        qqbar_randtest(x, state, 3, 50);
        qqbar_randtest(y, state, 3, 50);

        if (n_randint(state, 10) == 0)
        {
            qqbar_exp_pi_i(xr, n_randint(state, 10), 1 + n_randint(state, 3));
            qqbar_mul(y, x, xr);
        }

        qqbar_abs(xr, x);
        qqbar_abs(yr, y);

        c1 = qqbar_cmpabs(x, y);
        c2 = qqbar_cmp_re(xr, yr);

        if (c1 != c2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("xr = "); qqbar_print(xr); flint_printf("\n\n");
            flint_printf("yr = "); qqbar_print(yr); flint_printf("\n\n");
            flint_printf("%d\n\n", c1);
            flint_printf("%d\n\n", c2);
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(xr);
        qqbar_clear(yr);
    }

    TEST_FUNCTION_END(state);
}
