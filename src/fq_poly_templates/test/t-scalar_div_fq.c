/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_scalar_div_fq, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b;
        TEMPLATE(T, t) x;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, init) (x, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);
        while (TEMPLATE(T, is_zero(x, ctx)))
	    TEMPLATE(T, randtest(x, state, ctx));

        TEMPLATE(T, TEMPLATE(poly_scalar_div, T)) (b, a, x, ctx);
        TEMPLATE(T, TEMPLATE(poly_scalar_div, T)) (a, a, x, ctx);

        result = (TEMPLATE(T, poly_equal) (a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("x = "), TEMPLATE(T, print_pretty) (x, ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, clear) (x, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check a*x/x = a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b;
        TEMPLATE(T, t) x;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, init) (x, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);
        while (TEMPLATE(T, is_zero(x, ctx)))
	    TEMPLATE(T, randtest(x, state, ctx));

        TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (b, a, x, ctx);
        TEMPLATE(T, TEMPLATE(poly_scalar_div, T)) (b, b, x, ctx);

        result = (TEMPLATE(T, poly_equal) (a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("x = "), TEMPLATE(T, print_pretty) (x, ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, clear) (x, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
