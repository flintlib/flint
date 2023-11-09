/*
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

TEST_TEMPLATE_FUNCTION_START(T, poly_gcd, state)
{
    int i, result;

    /* Check that gcd(a,a) = a (made monic) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, g;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (g, ctx);

        TEMPLATE(T, poly_randtest_not_zero) (a, state, len, ctx);
        TEMPLATE(T, poly_make_monic) (b, a, ctx);
        TEMPLATE(T, poly_gcd) (g, a, a, ctx);

        result = (TEMPLATE(T, poly_equal) (g, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /*
       Find coprime polys, multiply by another poly
       and check the GCD is that poly
     */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong len, j;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c, g;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (g, ctx);

        TEMPLATE(T, poly_randtest_not_zero) (a, state, len, ctx);

        for (j = 0;
             (j < 100 * flint_test_multiplier())
             && !TEMPLATE(T, poly_is_one) (g, ctx); j++)
        {
            TEMPLATE(T, poly_randtest_not_zero) (b, state, len, ctx);
            TEMPLATE(T, poly_gcd) (g, a, b, ctx);

        }

        if (!TEMPLATE(T, poly_is_one) (g, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf
                ("could not find coprime polynomials after %wd tries\n",
                 j + 1);
            fflush(stdout);
            flint_abort();
        }

        for (j = 0; (j < 100 * flint_test_multiplier()) && (c->length < 2);
             j++)
            TEMPLATE(T, poly_randtest_not_zero) (c, state, len + 2, ctx);

        if (c->length < 2)
        {
            flint_printf("FAIL:\n");
            flint_printf
                ("could not find non-unit polynomial after %wd tries\n",
                 j + 1);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_make_monic) (c, c, ctx);

        TEMPLATE(T, poly_mul) (a, a, c, ctx);
        TEMPLATE(T, poly_mul) (b, b, c, ctx);

        TEMPLATE(T, poly_gcd) (g, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (g, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of a and g */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, g;

        len = n_randint(state, 200) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (g, ctx);

        TEMPLATE(T, poly_randtest_not_zero) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state, len, ctx);

        TEMPLATE(T, poly_gcd) (g, a, b, ctx);
        TEMPLATE(T, poly_gcd) (a, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (g, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of b and g */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, g;

        len = n_randint(state, 200) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (g, ctx);

        TEMPLATE(T, poly_randtest_not_zero) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state, len, ctx);

        TEMPLATE(T, poly_gcd) (g, a, b, ctx);
        TEMPLATE(T, poly_gcd) (b, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (g, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
