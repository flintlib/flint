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

TEST_FUNCTION_START(qqbar_abs, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, t, u;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(t);
        qqbar_init(u);

        qqbar_randtest(x, state, 3, 10);
        qqbar_randtest(y, state, 3, 10);

        qqbar_abs(z, x);
        qqbar_abs(t, y);
        qqbar_mul(t, t, z);

        qqbar_mul(u, x, y);
        qqbar_abs(u, u);

        if (!qqbar_equal(t, u) || !qqbar_is_real(t))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("t = "); qqbar_print(t); flint_printf("\n\n");
            flint_printf("u = "); qqbar_print(t); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(t);
        qqbar_clear(u);
    }

    TEST_FUNCTION_END(state);
}
