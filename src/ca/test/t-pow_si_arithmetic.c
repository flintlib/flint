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

TEST_FUNCTION_START(ca_pow_si_arithmetic, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, xa, xb, xaxb, xab;
        truth_t equal;
        slong a, b;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(xa, ctx);
        ca_init(xb, ctx);
        ca_init(xaxb, ctx);
        ca_init(xab, ctx);

        /* x^a * x^b = x^(a+b) */

        do {
            ca_randtest(x, state, 5, 5, ctx);
        } while (ca_check_is_zero(x, ctx) != T_FALSE);

        a = n_randint(state, 10) - 5;
        b = n_randint(state, 10) - 5;

        ca_pow_si_arithmetic(xa, x, a, ctx);
        ca_pow_si_arithmetic(xb, x, b, ctx);
        ca_mul(xaxb, xa, xb, ctx);
        ca_pow_si_arithmetic(xab, x, a + b, ctx);

        equal = ca_check_equal(xab, xaxb, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL x^a * x^b = x^(a+b)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("a = %wd", a); flint_printf("\n\n");
            flint_printf("b = %wd", b); flint_printf("\n\n");
            flint_printf("xa = "); ca_print(xa, ctx); flint_printf("\n\n");
            flint_printf("xb = "); ca_print(xb, ctx); flint_printf("\n\n");
            flint_printf("xaxb = "); ca_print(xaxb, ctx); flint_printf("\n\n");
            flint_printf("xab = "); ca_print(xab, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(xa, ctx);
        ca_clear(xb, ctx);
        ca_clear(xaxb, ctx);
        ca_clear(xab, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
