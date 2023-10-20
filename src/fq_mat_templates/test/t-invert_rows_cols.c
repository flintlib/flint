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

TEST_TEMPLATE_FUNCTION_START(T, mat_invert_rows_cols, state)
{
    slong m, n, rep, i, j;

    for (rep = 0; rep < 50 * flint_test_multiplier(); rep++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, mat_t) A;
        TEMPLATE(T, mat_t) B;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_init) (B, m, n, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);

        TEMPLATE(T, mat_set)(B, A, ctx);

        TEMPLATE(T, mat_invert_rows)(A, NULL, ctx);
        TEMPLATE(T, mat_invert_cols)(A, NULL, ctx);

        for (i = 0; i < A->r; i++)
        {
            for (j =0; j < A->c; j++)
            {
                if (!TEMPLATE(T, equal)(TEMPLATE(T, mat_entry)(B, i, j), TEMPLATE(T, mat_entry)(A, A->r - i - 1, A->c - j - 1), ctx))
                {
                    flint_printf("FAIL: B != A\n");
                    flint_printf("A:\n");
                    TEMPLATE(T, mat_print_pretty)(A, ctx);
                    flint_printf("B:\n");
                    TEMPLATE(T, mat_print_pretty)(B, ctx);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
