/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_equal_trunc, state)
{
    int i, result;

    /* Compare equal polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b;
        TEMPLATE(T, t) c;
        slong j, n;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 100), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 100), ctx);
        n = n_randint(state, 100);

        for (j = 0; j < n; j++)
        {
           TEMPLATE(T, poly_get_coeff) (c, a, j, ctx);
           TEMPLATE(T, poly_set_coeff) (b, j, c, ctx);
        }

        result = (TEMPLATE(T, poly_equal_trunc) (a, b, n, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear) (c, ctx);

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Compare unequal polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b;
        TEMPLATE(T, t) c, d;
        slong j, n, m;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, init) (c, ctx);
        TEMPLATE(T, init) (d, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 100), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 100), ctx);
        n = n_randint(state, 100) + 1;
        m = n_randint(state, n);

        for (j = 0; j < n; j++)
        {
           TEMPLATE(T, poly_get_coeff) (c, a, j, ctx);
           TEMPLATE(T, poly_set_coeff) (b, j, c, ctx);
        }

        TEMPLATE(T, poly_get_coeff) (c, a, m, ctx);
        TEMPLATE(T, set_si) (d, 1, ctx);
        TEMPLATE(T, add) (c, c, d, ctx);
        TEMPLATE(T, poly_set_coeff) (b, m, c, ctx);

        result = (!TEMPLATE(T, poly_equal_trunc) (a, b, n, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear) (c, ctx);
        TEMPLATE(T, clear) (d, ctx);

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
