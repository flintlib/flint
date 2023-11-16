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

TEST_FUNCTION_START(ca_sin_cos, state)
{
    slong iter;

    {
        ca_ctx_t ctx;
        ca_ctx_init(ctx);

        for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
        {
            ca_t x, s1, c1, s2, c2, s3, c3, s4, c4;

            if (n_randint(state, 10) == 0)
            {
                ca_ctx_clear(ctx);
                ca_ctx_init(ctx);
            }

            ca_init(x, ctx);
            ca_init(s1, ctx);
            ca_init(c1, ctx);
            ca_init(s2, ctx);
            ca_init(c2, ctx);
            ca_init(s3, ctx);
            ca_init(c3, ctx);
            ca_init(s4, ctx);
            ca_init(c4, ctx);

            ca_randtest_special(x, state, 5, 5, ctx);

            ca_sin_cos_direct(s1, c1, x, ctx);
            ca_sin_cos_exponential(s2, c2, x, ctx);
            ca_sin_cos_tangent(s3, c3, x, ctx);
            ca_sin_cos(s4, c4, x, ctx);

            if (ca_check_equal(s1, s2, ctx) == T_FALSE ||
                ca_check_equal(c1, c2, ctx) == T_FALSE ||
                ca_check_equal(s1, s3, ctx) == T_FALSE ||
                ca_check_equal(c1, c3, ctx) == T_FALSE ||
                ca_check_equal(s1, s4, ctx) == T_FALSE ||
                ca_check_equal(c1, c4, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
                flint_printf("s1 = "); ca_print(s1, ctx); flint_printf("\n\n");
                flint_printf("c1 = "); ca_print(c1, ctx); flint_printf("\n\n");
                flint_printf("s2 = "); ca_print(s2, ctx); flint_printf("\n\n");
                flint_printf("c2 = "); ca_print(c2, ctx); flint_printf("\n\n");
                flint_printf("s3 = "); ca_print(s3, ctx); flint_printf("\n\n");
                flint_printf("c3 = "); ca_print(c3, ctx); flint_printf("\n\n");
                flint_printf("s4 = "); ca_print(s4, ctx); flint_printf("\n\n");
                flint_printf("c4 = "); ca_print(c4, ctx); flint_printf("\n\n");
                flint_abort();
            }

            ca_clear(x, ctx);
            ca_clear(s1, ctx);
            ca_clear(c1, ctx);
            ca_clear(s2, ctx);
            ca_clear(c2, ctx);
            ca_clear(s3, ctx);
            ca_clear(c3, ctx);
            ca_clear(s4, ctx);
            ca_clear(c4, ctx);
        }

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
