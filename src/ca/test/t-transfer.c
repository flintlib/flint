/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca.h"

TEST_FUNCTION_START(ca_transfer, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx, ctx2;
        ca_t x, y, z;
        slong i, reps;

        ca_ctx_init(ctx);
        ca_ctx_init(ctx2);
        ca_init(x, ctx);
        ca_init(y, ctx2);
        ca_init(z, ctx);

        reps = 1 + n_randint(state, 10);
        for (i = 0; i < reps; i++)
        {
            ca_randtest_special(x, state, 5, 5, ctx);
            ca_randtest_special(y, state, 5, 5, ctx2);
            ca_randtest_special(z, state, 5, 5, ctx);

            ca_transfer(y, ctx2, x, ctx);
            ca_transfer(z, ctx, y, ctx2);

            if (ca_check_equal(x, z, ctx) == T_FALSE)
            {
                flint_printf("FAIL: not equal!\n");
                flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
                flint_printf("y = "); ca_print(y, ctx2); flint_printf("\n\n");
                flint_printf("z = "); ca_print(z, ctx); flint_printf("\n\n");
                flint_abort();
            }

            if (ca_check_equal(x, z, ctx) != T_TRUE && !ca_is_unknown(x, ctx))
            {
                flint_printf("FAIL: not equal!\n");
                flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
                flint_printf("y = "); ca_print(y, ctx2); flint_printf("\n\n");
                flint_printf("z = "); ca_print(z, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx2);
        ca_clear(z, ctx);
        ca_ctx_clear(ctx);
        ca_ctx_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
