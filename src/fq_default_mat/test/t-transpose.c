/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default_mat.h"

TEST_FUNCTION_START(fq_default_mat_transpose, state)
{
    slong m, n, rep;

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        fq_default_ctx_t ctx;
        fq_default_mat_t A, B, C;

        fq_default_ctx_init_randtest(ctx, state);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fq_default_mat_init(A, m, n, ctx);
        fq_default_mat_init(B, n, m, ctx);
        fq_default_mat_init(C, m, n, ctx);

        fq_default_mat_randtest(A, state, ctx);
        fq_default_mat_randtest(B, state, ctx);

        fq_default_mat_transpose(B, A, ctx);
        fq_default_mat_transpose(C, B, ctx);

        if (!fq_default_mat_equal(C, A, ctx))
        {
            flint_printf("FAIL: C != A\n");
            fflush(stdout);
            flint_abort();
        }

        fq_default_mat_clear(A, ctx);
        fq_default_mat_clear(B, ctx);
        fq_default_mat_clear(C, ctx);
        fq_default_ctx_clear(ctx);
    }

    /* Self-transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        fq_default_ctx_t ctx;
        fq_default_mat_t A, B;

        fq_default_ctx_init_randtest(ctx, state);

        m = n_randint(state, 20);

        fq_default_mat_init(A, m, m, ctx);
        fq_default_mat_init(B, m, m, ctx);

        fq_default_mat_randtest(A, state, ctx);
        fq_default_mat_set(B, A, ctx);
        fq_default_mat_transpose(B, B, ctx);
        fq_default_mat_transpose(B, B, ctx);

        if (!fq_default_mat_equal(B, A, ctx))
        {
            flint_printf("FAIL: B != A\n");
            fflush(stdout);
            flint_abort();
        }

        fq_default_mat_clear(A, ctx);
        fq_default_mat_clear(B, ctx);
        fq_default_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
