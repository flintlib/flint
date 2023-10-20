/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_solve, state)
{
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, mat_t) A, X, B, AX;
    slong i, m, n, r;
    int solved;

    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 50);
        n = n_randint(state, 50);

        TEMPLATE(T, mat_init) (A, m, m, ctx);
        TEMPLATE(T, mat_init) (B, m, n, ctx);
        TEMPLATE(T, mat_init) (X, m, n, ctx);
        TEMPLATE(T, mat_init) (AX, m, n, ctx);

        TEMPLATE(T, mat_randrank)(A, state, m, ctx);
        TEMPLATE(T, mat_randtest)(B, state, ctx);

        /* Dense */
        if (n_randint(state, 2))
            TEMPLATE(T, mat_randops)(A, 1+n_randint(state, 1+m*m), state, ctx);

        solved = TEMPLATE(T, mat_solve)(X, A, B, ctx);

        TEMPLATE(T, mat_mul)(AX, A, X, ctx);

        if (!TEMPLATE(T, mat_equal)(AX, B, ctx) || !solved)
        {
            flint_printf("FAIL:\n");
            flint_printf("AX != B!\n");
            flint_printf("A:\n");
            TEMPLATE(T, mat_print_pretty)(A, ctx);
            flint_printf("B:\n");
            TEMPLATE(T, mat_print_pretty)(B, ctx);
            flint_printf("X:\n");
            TEMPLATE(T, mat_print_pretty)(X, ctx);
            flint_printf("AX:\n");
            TEMPLATE(T, mat_print_pretty)(AX, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear)(A, ctx);
        TEMPLATE(T, mat_clear)(B, ctx);
        TEMPLATE(T, mat_clear)(X, ctx);
        TEMPLATE(T, mat_clear)(AX, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);
        r = n_randint(state, m);

        TEMPLATE(T, mat_init)(A, m, m, ctx);
        TEMPLATE(T, mat_init)(B, m, n, ctx);
        TEMPLATE(T, mat_init)(X, m, n, ctx);
        TEMPLATE(T, mat_init)(AX, m, n, ctx);

        TEMPLATE(T, mat_randrank)(A, state, r, ctx);
        TEMPLATE(T, mat_randtest)(B, state, ctx);

        /* Dense */
        if (n_randint(state, 2))
            TEMPLATE(T, mat_randops)(A, 1+n_randint(state, 1+m*m), state, ctx);

        solved = TEMPLATE(T, mat_solve)(X, A, B, ctx);

        if (solved)
        {
            flint_printf("FAIL:\n");
            flint_printf("singular system was 'solved'\n");
            TEMPLATE(T, mat_print_pretty)(A, ctx);
            TEMPLATE(T, mat_print_pretty)(X, ctx);
            TEMPLATE(T, mat_print_pretty)(B, ctx);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear)(A, ctx);
        TEMPLATE(T, mat_clear)(B, ctx);
        TEMPLATE(T, mat_clear)(X, ctx);
        TEMPLATE(T, mat_clear)(AX, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
