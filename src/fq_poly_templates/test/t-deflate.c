/*
    Copyright (C) 2011 Fredrik Johansson
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

TEST_TEMPLATE_FUNCTION_START(T, poly_deflate, state)
{
    int iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, poly_t) poly1, poly2, poly3;
        TEMPLATE(T, ctx_t) ctx;
        ulong infl1, infl, deflation;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (poly1, ctx);
        TEMPLATE(T, poly_init) (poly2, ctx);
        TEMPLATE(T, poly_init) (poly3, ctx);

        TEMPLATE(T, poly_randtest) (poly1, state, n_randint(state, 15), ctx);

        if (TEMPLATE(T, poly_length) (poly1, ctx) <= 1)
        {
            if (TEMPLATE(T, poly_deflation) (poly1, ctx) !=
                TEMPLATE(T, poly_length) (poly1, ctx))
            {
                flint_printf
                    ("FAIL: wrong deflation for constant polynomial\n");
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(T, poly_deflate) (poly2, poly1, n_randint(state, 5) + 1,
                                       ctx);
            if (!TEMPLATE(T, poly_equal) (poly2, poly1, ctx))
            {
                flint_printf
                    ("FAIL: constant polynomial changed on deflation\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {

            infl = n_randint(state, 13) + 1;
            infl1 = TEMPLATE(T, poly_deflation) (poly1, ctx);

            TEMPLATE(T, poly_inflate) (poly2, poly1, infl, ctx);

            deflation = TEMPLATE(T, poly_deflation) (poly2, ctx);

            if (deflation != infl * infl1)
            {
                flint_printf("FAIL: deflation = %wu, inflation: %wu, %wu\n",
                             deflation, infl, infl1);
                flint_printf("poly1:\n");
                TEMPLATE(T, poly_print) (poly1, ctx);
                flint_printf("\n\n");
                flint_printf("poly2:\n");
                TEMPLATE(T, poly_print) (poly2, ctx);
                flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(T, poly_deflate) (poly3, poly2, infl, ctx);
            if (!TEMPLATE(T, poly_equal) (poly3, poly1, ctx))
            {
                flint_printf("FAIL: deflation = %wu, inflation: %wu, %wu\n",
                             deflation, infl, infl1);
                flint_printf("Deflated polynomial not equal to input:\n");
                flint_printf("poly1:\n");
                TEMPLATE(T, poly_print) (poly1, ctx);
                flint_printf("\n\n");
                flint_printf("poly2:\n");
                TEMPLATE(T, poly_print) (poly2, ctx);
                flint_printf("\n\n");
                flint_printf("poly3:\n");
                TEMPLATE(T, poly_print) (poly3, ctx);
                flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(T, poly_deflate) (poly2, poly2, infl, ctx);
            if (!TEMPLATE(T, poly_equal) (poly3, poly2, ctx))
            {
                flint_printf("FAIL: aliasing\n");
                fflush(stdout);
                flint_abort();
            }
        }

        TEMPLATE(T, poly_clear) (poly1, ctx);
        TEMPLATE(T, poly_clear) (poly2, ctx);
        TEMPLATE(T, poly_clear) (poly3, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
