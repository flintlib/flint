/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca.h"

TEST_FUNCTION_START(ca_tan, state)
{
    slong iter;

    {
        ca_ctx_t ctx;
        ca_ctx_init(ctx);

        for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
        {
            ca_t x, t1, t2, t3, t4;

            if (n_randint(state, 10) == 0)
            {
                ca_ctx_clear(ctx);
                ca_ctx_init(ctx);
            }

            ca_init(x, ctx);
            ca_init(t1, ctx);
            ca_init(t2, ctx);
            ca_init(t3, ctx);
            ca_init(t4, ctx);

            ca_randtest_special(x, state, 5, 5, ctx);

            ca_tan_direct(t1, x, ctx);
            ca_tan_exponential(t2, x, ctx);
            ca_tan_sine_cosine(t3, x, ctx);
            ca_tan(t4, x, ctx);

            if (ca_check_equal(t1, t2, ctx) == T_FALSE ||
                ca_check_equal(t1, t3, ctx) == T_FALSE ||
                ca_check_equal(t1, t4, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
                flint_printf("t1 = "); ca_print(t1, ctx); flint_printf("\n\n");
                flint_printf("t2 = "); ca_print(t2, ctx); flint_printf("\n\n");
                flint_printf("t3 = "); ca_print(t3, ctx); flint_printf("\n\n");
                flint_printf("t4 = "); ca_print(t4, ctx); flint_printf("\n\n");
                flint_abort();
            }

            ca_clear(x, ctx);
            ca_clear(t1, ctx);
            ca_clear(t2, ctx);
            ca_clear(t3, ctx);
            ca_clear(t4, ctx);
        }

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
