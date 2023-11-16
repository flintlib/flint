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

TEST_FUNCTION_START(qqbar_root_ui, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, z;
        ulong n;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);

        qqbar_randtest(x, state, 4, 10);

        /* generate noise */
        qqbar_randtest(y, state, 2, 10);
        qqbar_add(x, x, y);
        qqbar_sub(x, x, y);

        n = 1 + n_randint(state, 4);

        qqbar_root_ui(y, x, n);

        if (n_randint(state, 2))
        {
            qqbar_pow_ui(z, y, n);
        }
        else
        {
            if (n == 1)
            {
                qqbar_set(z, y);
            }
            else if (n == 2)
            {
                qqbar_binary_op(z, y, y, 2);
            }
            else if (n == 3)
            {
                qqbar_binary_op(z, y, y, 2);
                qqbar_binary_op(z, z, y, 2);
            }
            else
            {
                qqbar_binary_op(z, y, y, 2);
                qqbar_binary_op(z, z, z, 2);
            }
        }

        if (!qqbar_equal(x, z) || (n > 1 && qqbar_csgn(y) == -1))
        {
            flint_printf("FAIL!\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
    }

    TEST_FUNCTION_END(state);
}
