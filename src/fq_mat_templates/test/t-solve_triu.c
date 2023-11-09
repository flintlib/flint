/*
    Copyright (C) 2010, 2011 Fredrik Johansson
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

TEST_TEMPLATE_FUNCTION_START(T, mat_solve_triu, state)
{
    slong i;

    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, X, B, Y;
        slong rows, cols;
        int unit;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        rows = n_randint(state, 50);
        cols = n_randint(state, 50);
        unit = n_randint(state, 2);

        TEMPLATE(T, mat_init) (A, rows, rows, ctx);
        TEMPLATE(T, mat_init) (B, rows, cols, ctx);
        TEMPLATE(T, mat_init) (X, rows, cols, ctx);
        TEMPLATE(T, mat_init) (Y, rows, cols, ctx);

        TEMPLATE(T, mat_randtriu) (A, state, unit, ctx);
        TEMPLATE(T, mat_randtest) (X, state, ctx);
        TEMPLATE(T, mat_mul) (B, A, X, ctx);

        /* Check Y = A^(-1) * (A * X) = X */
        TEMPLATE(T, mat_solve_triu) (Y, A, B, unit, ctx);
        if (!TEMPLATE(T, mat_equal) (Y, X, ctx))
        {
            printf("FAIL!\n");
            printf("A:\n");
            TEMPLATE(T, mat_print_pretty) (A, ctx);
            printf("X:\n");
            TEMPLATE(T, mat_print_pretty) (X, ctx);
            printf("B:\n");
            TEMPLATE(T, mat_print_pretty) (B, ctx);
            printf("Y:\n");
            TEMPLATE(T, mat_print_pretty) (Y, ctx);
            fflush(stdout);
            flint_abort();
        }

        /* Check aliasing */
        TEMPLATE(T, mat_solve_triu) (B, A, B, unit, ctx);
        if (!TEMPLATE(T, mat_equal) (B, X, ctx))
        {
            printf("FAIL!\n");
            printf("aliasing test failed");
            printf("A:\n");
            TEMPLATE(T, mat_print_pretty) (A, ctx);
            printf("B:\n");
            TEMPLATE(T, mat_print_pretty) (B, ctx);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (X, ctx);
        TEMPLATE(T, mat_clear) (Y, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
