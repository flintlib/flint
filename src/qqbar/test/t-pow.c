/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_pow, state)
{
    slong iter;

    /* Check x^m x^n = x^(m+n) */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, m, n, mn, xm, xn, xmxn, xmn;

        qqbar_init(x);
        qqbar_init(m);
        qqbar_init(n);
        qqbar_init(mn);
        qqbar_init(xm);
        qqbar_init(xn);
        qqbar_init(xmxn);
        qqbar_init(xmn);

        if (n_randint(state, 2))
        {
            qqbar_root_of_unity(x, n_randint(state, 100), 1 + n_randint(state, 10));
            qqbar_set_si(m, n_randtest(state));
            qqbar_div_ui(m, m, 1 + n_randint(state, 5));
            qqbar_set_si(n, n_randtest(state));
            qqbar_div_ui(n, n, 1 + n_randint(state, 5));
        }
        else
        {
            qqbar_randtest(x, state, 4, 10);
            qqbar_randtest(m, state, 4, 3);
            qqbar_randtest(n, state, 4, 3);
        }

        qqbar_add(mn, m, n);

        if (qqbar_pow(xm, x, m) && qqbar_pow(xn, x, n) && qqbar_pow(xmn, x, mn))
        {
            qqbar_mul(xmxn, xm, xn);

            if (!qqbar_equal(xmxn, xmn))
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("m = "); qqbar_print(m); flint_printf("\n\n");
                flint_printf("n = "); qqbar_print(n); flint_printf("\n\n");
                flint_printf("xm = "); qqbar_print(xm); flint_printf("\n\n");
                flint_printf("xn = "); qqbar_print(xn); flint_printf("\n\n");
                flint_printf("xmxn = "); qqbar_print(xmxn); flint_printf("\n\n");
                flint_printf("xmn = "); qqbar_print(xmn); flint_printf("\n\n");
                flint_abort();
            }
        }

        qqbar_clear(x);
        qqbar_clear(m);
        qqbar_clear(n);
        qqbar_clear(mn);
        qqbar_clear(xm);
        qqbar_clear(xn);
        qqbar_clear(xmxn);
        qqbar_clear(xmn);
    }

    TEST_FUNCTION_END(state);
}
