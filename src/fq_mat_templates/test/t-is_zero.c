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

TEST_TEMPLATE_FUNCTION_START(T, mat_is_zero, state)
{
    int i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A;
        TEMPLATE(T, t) x;
        slong j;

        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, mat_init) (A, rows, cols, ctx);

        if (!TEMPLATE(T, mat_is_zero) (A, ctx))
        {
            printf("FAIL!\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, init) (x, ctx);
        TEMPLATE(T, randtest_not_zero) (x, state, ctx);

        if (rows && cols)
        {
            j = n_randint(state, rows * cols);
            TEMPLATE(T, add) (A->entries + j, A->entries + j, x, ctx);

            if (TEMPLATE(T, mat_is_zero) (A, ctx))
            {
                printf("FAIL!\n");
                fflush(stdout);
                flint_abort();
            }
        }

        TEMPLATE(T, clear) (x, ctx);
        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
