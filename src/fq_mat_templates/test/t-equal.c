/*
    Copyright (C) 2010 William Hart
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

TEST_TEMPLATE_FUNCTION_START(T, mat_equal, state)
{
    int i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C, D, E;
        TEMPLATE(T, t) x;
        slong m, n, j;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, init) (x, ctx);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_init) (B, m, n, ctx);
        TEMPLATE(T, mat_init) (C, m, n, ctx);
        TEMPLATE(T, mat_init) (D, m + 1, n, ctx);
        TEMPLATE(T, mat_init) (E, m, n + 1, ctx);

        if (TEMPLATE(T, mat_equal) (A, D, ctx)
            || TEMPLATE(T, mat_equal) (A, E, ctx))
        {
            printf("FAIL: different dimensions should not be equal\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_set) (B, A, ctx);

        if (!TEMPLATE(T, mat_equal) (A, B, ctx))
        {
            printf("FAIL: copied matrices should be equal\n");
            fflush(stdout);
            flint_abort();
        }

        if (m && n)
        {
            j = n_randint(state, m * n);
            TEMPLATE(T, one) (x, ctx);
            TEMPLATE(T, add) (A->entries + j, A->entries + j, x, ctx);

            if (TEMPLATE(T, mat_equal) (A, B, ctx))
            {
                printf("FAIL: modified matrices should not be equal\n");
                fflush(stdout);
                flint_abort();
            }
        }

        TEMPLATE(T, clear) (x, ctx);
        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);
        TEMPLATE(T, mat_clear) (D, ctx);
        TEMPLATE(T, mat_clear) (E, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
