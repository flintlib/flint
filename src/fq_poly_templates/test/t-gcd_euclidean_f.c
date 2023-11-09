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

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_gcd_euclidean_f, state)
{
    int i, result;

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

    TEST_FUNCTION_END(state);
}
#endif
