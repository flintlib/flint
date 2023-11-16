/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_window_init_clear, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        slong r, c, r1, r2, c1, c2;
        gr_ctx_t ctx;
        gr_mat_t A, B, C;

        gr_ctx_init_random(ctx, state);

        r = n_randint(state, 8);
        c = n_randint(state, 8);

        r2 = n_randint(state, r + 1);
        r1 = n_randint(state, r2 + 1);

        c2 = n_randint(state, c + 1);
        c1 = n_randint(state, c2 + 1);

        gr_mat_init(A, r, c, ctx);
        gr_mat_window_init(B, A, r1, c1, r2, c2, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(B, state, ctx));

        gr_mat_window_init(C, A, r1, c1, r2, c2, ctx);

        status = gr_mat_one(C, ctx);

        if (status == GR_SUCCESS && gr_mat_is_one(C, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); gr_mat_print(C, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_mat_window_clear(B, ctx);
        gr_mat_window_clear(C, ctx);
        gr_mat_clear(A, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
