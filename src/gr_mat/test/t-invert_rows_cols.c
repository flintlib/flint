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

TEST_FUNCTION_START(gr_mat_invert_rows_cols, state)
{
    slong iter;

    for (iter = 0; iter < 100; iter++)
    {
        slong m, n, i, j;
        gr_ctx_t ctx;
        gr_mat_t A;
        gr_mat_t B;

        gr_ctx_init_random(ctx, state);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        gr_mat_init(A, m, n, ctx);
        gr_mat_init(B, m, n, ctx);

        if (gr_mat_randtest(A, state, ctx) == GR_SUCCESS &&
            gr_mat_set(B, A, ctx) == GR_SUCCESS)
        {
            GR_MUST_SUCCEED(gr_mat_invert_rows(A, NULL, ctx));
            GR_MUST_SUCCEED(gr_mat_invert_cols(A, NULL, ctx));

            for (i = 0; i < A->r; i++)
            {
                for (j = 0; j < A->c; j++)
                {
                    if (gr_equal(gr_mat_entry_ptr(B, i, j, ctx), gr_mat_entry_ptr(A, A->r - i - 1, A->c - j - 1, ctx), ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: B != A\n");
                        flint_printf("A:\n");
                        gr_mat_print(A, ctx);
                        flint_printf("B:\n");
                        gr_mat_print(B, ctx);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
