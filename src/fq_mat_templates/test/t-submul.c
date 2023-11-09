/*
    Copyright (C) 2010 Fredrik Johansson
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

TEST_TEMPLATE_FUNCTION_START(T, mat_submul, state)
{
    slong i;

    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C, D, T, E;

        slong m, k, n;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 50);
        k = n_randint(state, 50);
        n = n_randint(state, 50);

        TEMPLATE(T, mat_init) (A, m, k, ctx);
        TEMPLATE(T, mat_init) (B, k, n, ctx);
        TEMPLATE(T, mat_init) (C, m, n, ctx);
        TEMPLATE(T, mat_init) (D, m, n, ctx);
        TEMPLATE(T, mat_init) (T, m, n, ctx);
        TEMPLATE(T, mat_init) (E, m, n, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_randtest) (B, state, ctx);
        TEMPLATE(T, mat_randtest) (C, state, ctx);

        TEMPLATE(T, mat_submul) (D, C, A, B, ctx);

        TEMPLATE(T, mat_mul) (T, A, B, ctx);
        TEMPLATE(T, mat_sub) (E, C, T, ctx);

        if (!TEMPLATE(T, mat_equal) (D, E, ctx))
        {
            printf("FAIL: results not equal\n");
            TEMPLATE(T, mat_print_pretty) (A, ctx);
            TEMPLATE(T, mat_print_pretty) (B, ctx);
            TEMPLATE(T, mat_print_pretty) (C, ctx);
            TEMPLATE(T, mat_print_pretty) (D, ctx);
            TEMPLATE(T, mat_print_pretty) (E, ctx);
            fflush(stdout);
            flint_abort();
        }

        /* Check aliasing */
        TEMPLATE(T, mat_submul) (C, C, A, B, ctx);

        if (!TEMPLATE(T, mat_equal) (C, E, ctx))
        {
            printf("FAIL: results not equal (aliasing)\n");
            TEMPLATE(T, mat_print_pretty) (A, ctx);
            TEMPLATE(T, mat_print_pretty) (B, ctx);
            TEMPLATE(T, mat_print_pretty) (C, ctx);
            TEMPLATE(T, mat_print_pretty) (D, ctx);
            TEMPLATE(T, mat_print_pretty) (E, ctx);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);
        TEMPLATE(T, mat_clear) (D, ctx);
        TEMPLATE(T, mat_clear) (E, ctx);
        TEMPLATE(T, mat_clear) (T, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
