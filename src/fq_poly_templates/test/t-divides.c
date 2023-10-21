/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_divides, state)
{
    int i, result;

    /* Check that b divides a b and that the quotient is a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c, q;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (q, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 50) + 1, ctx);
        TEMPLATE(T, poly_mul) (c, a, b, ctx);

        result = TEMPLATE(T, poly_divides) (q, c, b, ctx)
            && TEMPLATE(T, poly_equal) (q, a, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            flint_printf("q = "), TEMPLATE(T, poly_print_pretty) (q, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of a and q */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 50) + 1, ctx);
        TEMPLATE(T, poly_mul) (c, a, b, ctx);

        result = TEMPLATE(T, poly_divides) (c, c, b, ctx)
            && TEMPLATE(T, poly_equal) (c, a, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of b and q */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 50) + 1, ctx);
        TEMPLATE(T, poly_mul) (c, a, b, ctx);

        result = TEMPLATE(T, poly_divides) (b, c, b, ctx)
            && TEMPLATE(T, poly_equal) (b, a, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
