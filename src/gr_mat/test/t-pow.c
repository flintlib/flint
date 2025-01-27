/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_pow, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        slong n, i;
        ulong e;
        gr_ctx_t ctx;
        gr_mat_t A, B, C;

        gr_ctx_init_random(ctx, state);

        n = n_randint(state, 6);
        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C, n, n, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(B, state, ctx));

        GR_MUST_SUCCEED(gr_mat_pow_ui(B, A, e, ctx));

        GR_MUST_SUCCEED(gr_mat_one(C, ctx));
        for (i = 0; i < e; i++)
            GR_MUST_SUCCEED(gr_mat_mul(C, C, A, ctx));

        if (!gr_mat_equal(C, B, ctx))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_mat_pow_ui(A, A, e, ctx));

        if (!gr_mat_equal(A, B, ctx))
        {
            flint_printf("FAIL: aliasing failed\n");
            fflush(stdout);
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
