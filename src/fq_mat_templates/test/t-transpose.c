/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_transpose, state)
{
    slong m, n, rep;

    /* Rectangular transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        TEMPLATE(T, mat_t) A, B, C;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_init_randtest) (ctx, state, 3);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_init) (B, n, m, ctx);
        TEMPLATE(T, mat_init) (C, m, n, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_randtest) (B, state, ctx);

        TEMPLATE(T, mat_transpose) (B, A, ctx);
        TEMPLATE(T, mat_transpose) (C, B, ctx);

        if (!TEMPLATE(T, mat_equal) (C, A, ctx))
        {
            flint_printf("FAIL: C != A\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Self-transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        TEMPLATE(T, mat_t) A, B;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_init_randtest) (ctx, state, 3);

        m = n_randint(state, 20);

        TEMPLATE(T, mat_init) (A, m, m, ctx);
        TEMPLATE(T, mat_init) (B, m, m, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_set) (B, A, ctx);
        TEMPLATE(T, mat_transpose) (B, B, ctx);
        TEMPLATE(T, mat_transpose) (B, B, ctx);

        if (!TEMPLATE(T, mat_equal) (B, A, ctx))
        {
            flint_printf("FAIL: B != A\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
