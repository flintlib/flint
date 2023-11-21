/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_rising_ui, state)
{
    slong iter;

    /* check functional equation */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, xn, y, z;
        ulong n, m;

        acb_init(x);
        acb_init(xn);
        acb_init(y);
        acb_init(z);

        acb_randtest(x, state, 1 + n_randint(state, 4000), 10);
        n = n_randint(state, 80);
        m = n_randint(state, 40);
        acb_add_ui(xn, x, n, 1 + n_randint(state, 4000));

        acb_rising_ui(y, x, n, 2 + n_randint(state, 4000));
        acb_rising_ui(z, xn, m, 2 + n_randint(state, 4000));
        acb_mul(y, y, z, 2 + n_randint(state, 4000));

        acb_rising_ui(z, x, n + m, 2 + n_randint(state, 4000));

        if (!acb_overlaps(y, z))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu\n", n);
            flint_printf("m = %wu\n", m);
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("xn = "); acb_print(xn); flint_printf("\n\n");
            flint_printf("y = "); acb_print(y); flint_printf("\n\n");
            flint_printf("z = "); acb_print(z); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(xn);
        acb_clear(y);
        acb_clear(z);
    }

    /* aliasing of y and x */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, y;
        ulong n;
        slong prec;

        acb_init(x);
        acb_init(y);

        acb_randtest(x, state, 1 + n_randint(state, 200), 10);
        acb_randtest(y, state, 1 + n_randint(state, 200), 10);
        n = n_randint(state, 100);

        prec = 2 + n_randint(state, 4000);
        acb_rising_ui(y, x, n, prec);
        acb_rising_ui(x, x, n, prec);

        if (!acb_equal(x, y))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("y = "); acb_print(y); flint_printf("\n\n");
            flint_printf("n = %wu\n", n);
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y);
    }

    TEST_FUNCTION_END(state);
}
