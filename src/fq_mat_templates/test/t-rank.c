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

TEST_TEMPLATE_FUNCTION_START(T, mat_rank, state)
{
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, mat_t) A;
    slong i, m, n, d, r;

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        TEMPLATE(T, ctx_randtest) (ctx, state);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            TEMPLATE(T, mat_init) (A, m, n, ctx);
            TEMPLATE(T, mat_randrank) (A, state, r, ctx);
            if (r != TEMPLATE(T, mat_rank) (A, ctx))
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }
            TEMPLATE(T, mat_clear) (A, ctx);
        }

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Dense */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        TEMPLATE(T, ctx_randtest) (ctx, state);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            d = n_randint(state, 2 * m * n + 1);
            TEMPLATE(T, mat_init) (A, m, n, ctx);
            TEMPLATE(T, mat_randrank) (A, state, r, ctx);
            TEMPLATE(T, mat_randops) (A, d, state, ctx);
            if (r != TEMPLATE(T, mat_rank) (A, ctx))
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }
            TEMPLATE(T, mat_clear) (A, ctx);
        }
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
