/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart
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

    flint_printf("divrem_divconquer....");
    fflush(stdout);

    /* Check q*b + r = a, no aliasing */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, q, r, t;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_init) (t, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 100), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 100) + 1, ctx);

        TEMPLATE(T, poly_divrem_divconquer) (q, r, a, b, ctx);

        TEMPLATE(T, poly_mul) (t, q, b, ctx);
        TEMPLATE(T, poly_add) (t, t, r, ctx);

        result = (TEMPLATE(T, poly_equal) (a, t, ctx));
        if (!result)
        {
            flint_printf("FAIL #1:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print) (a, ctx),
                flint_printf("\n\n");
            flint_printf("b = "), TEMPLATE(T, poly_print) (b, ctx),
                flint_printf("\n\n");
            flint_printf("q = "), TEMPLATE(T, poly_print) (q, ctx),
                flint_printf("\n\n");
            flint_printf("r = "), TEMPLATE(T, poly_print) (r, ctx),
                flint_printf("\n\n");
            flint_printf("t = "), TEMPLATE(T, poly_print) (t, ctx),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Alias a and q, b and r */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, q, r;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 100), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 100) + 1, ctx);

        TEMPLATE(T, poly_divrem_divconquer) (q, r, a, b, ctx);
        TEMPLATE(T, poly_divrem_divconquer) (a, b, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (q, a, ctx)
                  && TEMPLATE(T, poly_equal) (r, b, ctx));
        if (!result)
        {
            flint_printf("FAIL #2:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print) (a, ctx),
                flint_printf("\n\n");
            flint_printf("b = "), TEMPLATE(T, poly_print) (b, ctx),
                flint_printf("\n\n");
            flint_printf("q = "), TEMPLATE(T, poly_print) (q, ctx),
                flint_printf("\n\n");
            flint_printf("r = "), TEMPLATE(T, poly_print) (r, ctx),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Alias b and q, a and r */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, q, r;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 100), ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 100) + 1, ctx);

        TEMPLATE(T, poly_divrem_divconquer) (q, r, a, b, ctx);
        TEMPLATE(T, poly_divrem_divconquer) (b, a, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (q, b, ctx)
                  && TEMPLATE(T, poly_equal) (r, a, ctx));
        if (!result)
        {
            flint_printf("FAIL #3:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print) (a, ctx),
                flint_printf("\n\n");
            flint_printf("b = "), TEMPLATE(T, poly_print) (b, ctx),
                flint_printf("\n\n");
            flint_printf("q = "), TEMPLATE(T, poly_print) (q, ctx),
                flint_printf("\n\n");
            flint_printf("r = "), TEMPLATE(T, poly_print) (r, ctx),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
