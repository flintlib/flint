/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
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

TEST_TEMPLATE_FUNCTION_START(T, poly_sub, state)
{
    int i, result;

    /* Check aliasing: a = a - b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest) (b, state, len, ctx);

        TEMPLATE(T, poly_sub) (c, a, b, ctx);
        TEMPLATE(T, poly_sub) (a, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (aliasing a = a - b):\n\n");
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

    /* Check aliasing: b = a - b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest) (b, state, len, ctx);

        TEMPLATE(T, poly_sub) (c, a, b, ctx);
        TEMPLATE(T, poly_sub) (b, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (b = a - b):\n\n");
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

    /* Check aliasing: a = a - a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, c;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);

        TEMPLATE(T, poly_sub) (c, a, a, ctx);
        TEMPLATE(T, poly_sub) (a, a, a, ctx);

        result = (TEMPLATE(T, poly_equal) (a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (a = a - a):\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check that a - b == -(b - a) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c, e;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (e, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest) (b, state, len, ctx);

        TEMPLATE(T, poly_sub) (c, a, b, ctx);
        TEMPLATE(T, poly_sub) (e, b, a, ctx);
        TEMPLATE(T, poly_neg) (e, e, ctx);

        result = (TEMPLATE(T, poly_equal) (e, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (a - b == -(b - a)):\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            flint_printf("e = "), TEMPLATE(T, poly_print_pretty) (e, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (e, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
