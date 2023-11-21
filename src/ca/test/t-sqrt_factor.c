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

TEST_FUNCTION_START(ca_sqrt_factor, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, z, a, b, c, d;
        truth_t equal;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(z, ctx);
        ca_init(a, ctx);
        ca_init(b, ctx);
        ca_init(c, ctx);
        ca_init(d, ctx);

        /* test sqrt(x)^2 = x */
        ca_randtest_special(x, state, 5, 5, ctx);
        ca_randtest_special(y, state, 5, 5, ctx);
        ca_randtest_special(z, state, 5, 5, ctx);
        ca_randtest_special(a, state, 5, 5, ctx);
        ca_randtest_special(b, state, 5, 5, ctx);

        ca_mul(x, x, x, ctx);
        ca_mul(x, x, z, ctx);

        /* todo: random flags */
        ca_sqrt_factor(y, x, CA_FACTOR_ZZ_SMOOTH | CA_FACTOR_POLY_FULL, ctx);
        ca_mul(z, y, y, ctx);
        ca_sub(b, x, z, ctx);

        equal = ca_check_equal(x, z, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL (sqrt(x)^2 != x)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_printf("z = "); ca_print(z, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_sqrt_inert(z, x, ctx);

        equal = ca_check_equal(y, z, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL: sqrt(x) != sqrt_inert(x)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_printf("z = "); ca_print(z, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(z, ctx);
        ca_clear(a, ctx);
        ca_clear(b, ctx);
        ca_clear(c, ctx);
        ca_clear(d, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
