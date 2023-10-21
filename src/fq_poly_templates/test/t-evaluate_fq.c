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

TEST_TEMPLATE_FUNCTION_START(T, poly_evaluate_fq, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) f;
        TEMPLATE(T, t) x, y, z;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, init) (x, ctx);
        TEMPLATE(T, init) (y, ctx);
        TEMPLATE(T, init) (z, ctx);

        TEMPLATE(T, poly_randtest) (f, state, len, ctx);
        TEMPLATE(T, randtest) (x, state, ctx);

        TEMPLATE(T, set) (z, x, ctx);
        TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (y, f, x, ctx);
        TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (x, f, x, ctx);

        result = (TEMPLATE(T, equal) (x, y, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("f = "), TEMPLATE(T, poly_print_pretty) (f, "X", ctx),
                flint_printf("\n");
            flint_printf("x = "), TEMPLATE(T, print_pretty) (x, ctx),
                flint_printf("\n");
            flint_printf("y = "), TEMPLATE(T, print_pretty) (y, ctx),
                flint_printf("\n");
            flint_printf("z = "), TEMPLATE(T, print_pretty) (z, ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, clear) (x, ctx);
        TEMPLATE(T, clear) (y, ctx);
        TEMPLATE(T, clear) (z, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check (f + g)(a) == f(a) + g(a) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) f, g, h;
        TEMPLATE(T, t) x, y, z;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (h, ctx);
        TEMPLATE(T, init) (x, ctx);
        TEMPLATE(T, init) (y, ctx);
        TEMPLATE(T, init) (z, ctx);

        TEMPLATE(T, poly_randtest) (f, state, len, ctx);
        TEMPLATE(T, poly_randtest) (g, state, len, ctx);
        TEMPLATE(T, randtest) (x, state, ctx);

        TEMPLATE(T, poly_add) (h, f, g, ctx);
        TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (y, f, x, ctx);
        TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (z, g, x, ctx);
        TEMPLATE(T, add) (y, y, z, ctx);
        TEMPLATE(T, TEMPLATE(poly_evaluate, T)) (z, h, x, ctx);

        result = (TEMPLATE(T, equal) (y, z, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("f = "), TEMPLATE(T, poly_print_pretty) (f, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
            flint_printf("h = "), TEMPLATE(T, poly_print_pretty) (h, "X", ctx),
                flint_printf("\n");
            flint_printf("x = "), TEMPLATE(T, print_pretty) (x, ctx),
                flint_printf("\n");
            flint_printf("y = "), TEMPLATE(T, print_pretty) (y, ctx),
                flint_printf("\n");
            flint_printf("z = "), TEMPLATE(T, print_pretty) (z, ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (h, ctx);
        TEMPLATE(T, clear) (x, ctx);
        TEMPLATE(T, clear) (y, ctx);
        TEMPLATE(T, clear) (z, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
