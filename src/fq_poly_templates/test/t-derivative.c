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

TEST_TEMPLATE_FUNCTION_START(T, poly_derivative, state)
{
    int i, result;

    /* Check aliasing */
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
        TEMPLATE(T, poly_set) (b, a, ctx);
        TEMPLATE(T, poly_derivative) (c, b, ctx);
        TEMPLATE(T, poly_derivative) (b, b, ctx);

        result = (TEMPLATE(T, poly_equal) (b, c, ctx));
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

    /* Check constants have derivative zero */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 2), ctx);
        TEMPLATE(T, poly_derivative) (b, a, ctx);

        result = (TEMPLATE(T, poly_is_zero) (b, ctx));
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

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check (f g)' == f' g + f g'  */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c, d, lhs, rhs;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_init) (lhs, ctx);
        TEMPLATE(T, poly_init) (rhs, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 100), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 100), ctx);

        TEMPLATE(T, poly_mul) (lhs, a, b, ctx);
        TEMPLATE(T, poly_derivative) (lhs, lhs, ctx);
        TEMPLATE(T, poly_derivative) (c, a, ctx);
        TEMPLATE(T, poly_derivative) (d, b, ctx);
        TEMPLATE(T, poly_mul) (c, c, b, ctx);
        TEMPLATE(T, poly_mul) (d, a, d, ctx);
        TEMPLATE(T, poly_add) (rhs, c, d, ctx);

        result = (TEMPLATE(T, poly_equal) (lhs, rhs, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a   = "), TEMPLATE(T, poly_print_pretty) (a, "X",
                                                                    ctx),
                flint_printf("\n");
            flint_printf("b   = "), TEMPLATE(T, poly_print_pretty) (b, "X",
                                                                    ctx),
                flint_printf("\n");
            flint_printf("c   = "), TEMPLATE(T, poly_print_pretty) (c, "X",
                                                                    ctx),
                flint_printf("\n");
            flint_printf("d   = "), TEMPLATE(T, poly_print_pretty) (d, "X",
                                                                    ctx),
                flint_printf("\n");
            flint_printf("lhs = "), TEMPLATE(T, poly_print_pretty) (lhs, "X",
                                                                    ctx),
                flint_printf("\n");
            flint_printf("rhs = "), TEMPLATE(T, poly_print_pretty) (rhs, "X",
                                                                    ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, poly_clear) (lhs, ctx);
        TEMPLATE(T, poly_clear) (rhs, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
