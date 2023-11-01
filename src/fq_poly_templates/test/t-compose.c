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

TEST_TEMPLATE_FUNCTION_START(T, poly_compose, state)
{
    int i, result;

    /* Check aliasing of the first argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) f, g, h;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (h, ctx);

        TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 40), ctx);
        TEMPLATE(T, poly_randtest) (g, state, n_randint(state, 20), ctx);

        TEMPLATE(T, poly_compose) (h, f, g, ctx);
        TEMPLATE(T, poly_compose) (f, f, g, ctx);

        result = (TEMPLATE(T, poly_equal) (f, h, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("f = "), TEMPLATE(T, poly_print_pretty) (f, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
            flint_printf("h = "), TEMPLATE(T, poly_print_pretty) (h, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (h, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) f, g, h;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (h, ctx);

        TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 40), ctx);
        TEMPLATE(T, poly_randtest) (g, state, n_randint(state, 20), ctx);

        TEMPLATE(T, poly_compose) (h, f, g, ctx);
        TEMPLATE(T, poly_compose) (g, f, g, ctx);

        result = (TEMPLATE(T, poly_equal) (g, h, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("f = "), TEMPLATE(T, poly_print_pretty) (f, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
            flint_printf("h = "), TEMPLATE(T, poly_print_pretty) (h, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (h, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Compare with the naive method */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) f, g, h, s, t;
        slong k;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (h, ctx);
        TEMPLATE(T, poly_init) (s, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (g, state, n_randint(state, 40), ctx);
        TEMPLATE(T, poly_randtest) (h, state, n_randint(state, 20), ctx);

        TEMPLATE(T, poly_one) (t, ctx);
        for (k = 0; k < TEMPLATE(T, poly_length) (g, ctx); k++)
        {
            TEMPLATE(T, TEMPLATE(poly_scalar_addmul, T)) (s, t, g->coeffs + k,
                                                          ctx);
            TEMPLATE(T, poly_mul) (t, t, h, ctx);
        }

        TEMPLATE(T, poly_compose) (f, g, h, ctx);

        result = (TEMPLATE(T, poly_equal) (f, s, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("f = "), TEMPLATE(T, poly_print_pretty) (f, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
            flint_printf("h = "), TEMPLATE(T, poly_print_pretty) (h, "X", ctx),
                flint_printf("\n");
            flint_printf("s = "), TEMPLATE(T, poly_print_pretty) (s, "X", ctx),
                flint_printf("\n");
            flint_printf("t = "), TEMPLATE(T, poly_print_pretty) (t, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (h, ctx);
        TEMPLATE(T, poly_clear) (s, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
