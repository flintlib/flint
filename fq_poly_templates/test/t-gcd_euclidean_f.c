/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_euclidean_f....");
    fflush(stdout);

    /*
       Compare with the usual GCD function.

       N.B.  I checked by hand that this test shows both outcomes, 
       i.e. trivial and non-trivial factors, sufficiently frequently.
     */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) f;
        TEMPLATE(T, poly_t) a, b, c, d;

        TEMPLATE(T, ctx_randtest_reducible) (ctx, state);
        TEMPLATE(T, init) (f, ctx);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 60), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 60), ctx);

        TEMPLATE(T, poly_gcd_euclidean_f) (f, c, a, b, ctx);
        if (!TEMPLATE(T, is_one) (f, ctx))
        {
            result = 1;
        }
        else
        {
            TEMPLATE(T, poly_gcd) (d, a, b, ctx);
            result = TEMPLATE(T, poly_equal) (c, d, ctx);
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "x", ctx);
            flint_printf("\n\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "x", ctx);
            flint_printf("\n\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "x", ctx);
            flint_printf("\n\n");
            flint_printf("d = "), TEMPLATE(T, poly_print_pretty) (d, "x", ctx);
            flint_printf("\n\n");
            flint_printf("f = "), TEMPLATE(T, print_pretty) (f, ctx),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, clear) (f, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);

    /* Check aliasing of a and g */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, g;
        TEMPLATE(T, t) f1, f2;

        len = n_randint(state, 60) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, init) (f1, ctx);
        TEMPLATE(T, init) (f2, ctx);

        TEMPLATE(T, poly_randtest_not_zero) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state, len, ctx);

        TEMPLATE(T, poly_gcd_euclidean_f) (f1, g, a, b, ctx);
        TEMPLATE(T, poly_gcd_euclidean_f) (f2, a, a, b, ctx);

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
            flint_printf("f1 = "), TEMPLATE(T, print) (f1, ctx), flint_printf("\n");
	    flint_printf("f2 = "), TEMPLATE(T, print) (f2, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, clear) (f1, ctx);
        TEMPLATE(T, clear) (f2, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }


    /* Check aliasing of b and g */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, g;
        TEMPLATE(T, t) f1, f2;

        len = n_randint(state, 60) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, init) (f1, ctx);
        TEMPLATE(T, init) (f2, ctx);

        TEMPLATE(T, poly_randtest_not_zero) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state, len, ctx);

        TEMPLATE(T, poly_gcd_euclidean_f) (f1, g, a, b, ctx);
        TEMPLATE(T, poly_gcd_euclidean_f) (f2, b, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (g, b, ctx) && TEMPLATE(T, equal) (f1, f2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf(":wq\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("g = "), TEMPLATE(T, poly_print_pretty) (g, "X", ctx),
                flint_printf("\n");
	    flint_printf("f1 = "), TEMPLATE(T, print) (f1, ctx), flint_printf("\n");
            flint_printf("f2 = "), TEMPLATE(T, print) (f2, ctx), flint_printf("\n");
	    fflush(stdout);
	    flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, clear) (f1, ctx);
        TEMPLATE(T, clear) (f2, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    flint_printf("PASS\n");
    return 0;
}



#endif
