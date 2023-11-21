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

TEST_FUNCTION_START(ca_exp, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, z, a, b, c, d, e;
        truth_t equal;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(z, ctx);
        ca_init(a, ctx);
        ca_init(b, ctx);
        ca_init(c, ctx);
        ca_init(d, ctx);
        ca_init(e, ctx);

        /* exp(x+y+z) = exp(x)*exp(y)*exp(z) */

        ca_randtest(x, state, 5, 5, ctx);
        ca_randtest(y, state, 5, 5, ctx);
        ca_randtest(z, state, 5, 5, ctx);

        if (n_randint(state, 4) == 0)
        {
            if (n_randint(state, 2))
            {
                do {
                    ca_randtest(a, state, 5, 5, ctx);
                    ca_log(a, a, ctx);
                } while (ca_check_is_number(a, ctx) != T_TRUE);
                ca_mul(x, x, a, ctx);
            }
            else
            {
                ca_pi_i(a, ctx);
                ca_mul(x, x, a, ctx);
            }
        }

        ca_exp(a, x, ctx);
        ca_exp(b, y, ctx);
        ca_exp(c, z, ctx);

        ca_add(d, x, y, ctx);
        ca_add(d, d, z, ctx);
        ca_exp(d, d, ctx);

        ca_mul(e, a, b, ctx);
        ca_mul(e, e, c, ctx);

        equal = ca_check_equal(d, e, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL exp(x+y+z) = exp(x)*exp(y)*exp(z)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_printf("z = "); ca_print(z, ctx); flint_printf("\n\n");
            flint_printf("a = "); ca_print(a, ctx); flint_printf("\n\n");
            flint_printf("b = "); ca_print(b, ctx); flint_printf("\n\n");
            flint_printf("c = "); ca_print(c, ctx); flint_printf("\n\n");
            flint_printf("d = "); ca_print(d, ctx); flint_printf("\n\n");
            flint_printf("e = "); ca_print(e, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(z, ctx);
        ca_clear(a, ctx);
        ca_clear(b, ctx);
        ca_clear(c, ctx);
        ca_clear(d, ctx);
        ca_clear(e, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
