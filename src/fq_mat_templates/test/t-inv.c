/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_inv, state)
{
    TEMPLATE(T, mat_t) A, B, C, I;
    TEMPLATE(T, ctx_t) ctx;
    slong i, j, m, r;
    int result;

    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, mat_init)(A, m, m, ctx);
        TEMPLATE(T, mat_init)(B, m, m, ctx);
        TEMPLATE(T, mat_init)(C, m, m, ctx);
        TEMPLATE(T, mat_init)(I, m, m, ctx);

        for (j = 0; j < m; j++)
            TEMPLATE(T, one)(TEMPLATE(T, mat_entry)(I, j, j), ctx);

        /* Verify that A * A^-1 = I for random matrices */

        TEMPLATE(T, mat_randrank)(A, state, m, ctx);
        /* Dense or sparse? */
        if (n_randint(state, 2))
            TEMPLATE(T, mat_randops)(A, 1+n_randint(state, 1+m*m), state, ctx);

        result = TEMPLATE(T, mat_inv)(B, A, ctx);
        TEMPLATE(T, mat_mul)(C, A, B, ctx);

        if (!TEMPLATE(T, mat_equal)(C, I, ctx) || !result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A * A^-1 != I!\n");
            flint_printf("A:\n");
            TEMPLATE(T, mat_print_pretty)(A, ctx);
            flint_printf("A^-1:\n");
            TEMPLATE(T, mat_print_pretty)(B, ctx);
            flint_printf("A * A^-1:\n");
            TEMPLATE(T, mat_print_pretty)(C, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* Test aliasing */
        TEMPLATE(T, mat_set)(C, A, ctx);
        TEMPLATE(T, mat_inv)(A, A, ctx);
        TEMPLATE(T, mat_mul)(B, A, C, ctx);

        if (!TEMPLATE(T, mat_equal)(B, I, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf("aliasing failed!\n");
            TEMPLATE(T, mat_print_pretty)(C, ctx);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear)(A, ctx);
        TEMPLATE(T, mat_clear)(B, ctx);
        TEMPLATE(T, mat_clear)(C, ctx);
        TEMPLATE(T, mat_clear)(I, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        r = n_randint(state, m);
        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, mat_init)(A, m, m, ctx);
        TEMPLATE(T, mat_init)(B, m, m, ctx);

        TEMPLATE(T, mat_randrank)(A, state, r, ctx);

        /* Dense */
        if (n_randint(state, 2))
            TEMPLATE(T, mat_randops)(A, 1+n_randint(state, 1+m*m), state, ctx);

        result = TEMPLATE(T, mat_inv)(B, A, ctx);

        if (result)
        {
            flint_printf("FAIL:\n");
            flint_printf("singular matrix reported as invertible\n");
            fflush(stdout);
            flint_abort();
        }

        /* Aliasing */
        result = TEMPLATE(T, mat_inv)(A, A, ctx);
        if (result)
        {
            flint_printf("FAIL:\n");
            flint_printf("singular matrix reported as invertible\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear)(A, ctx);
        TEMPLATE(T, mat_clear)(B, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
