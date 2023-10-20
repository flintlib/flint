/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_one, state)
{
    int iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A;
        slong m, n;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_one) (A, ctx);

        if (!TEMPLATE(T, mat_is_one) (A, ctx))
        {
            printf("FAIL: expected matrix to be one\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
