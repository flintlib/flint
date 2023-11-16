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

TEST_GR_FUNCTION_START(gr_mat_nullspace, state, count_success, count_unable, count_domain)
{
    slong iter;

    /* Check that random row/column operations preserve rank */
    for (iter = 0; iter < 10000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, X, AX;
        slong r, c, rank, nullity;
        int status = GR_SUCCESS;
        int status2 = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        gr_mat_init(A, r, c, ctx);
        gr_mat_init(X, 0, 0, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

        status = gr_mat_nullspace(X, A, ctx);

        if (status == GR_SUCCESS)
        {
            status2 = gr_mat_rank(&rank, A, ctx);

            if (status2 == GR_SUCCESS)
            {
                nullity = c - rank;

                if (nullity != gr_mat_ncols(X, ctx))
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("r = %wd, c = %wd\n", r, c);
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); gr_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, nullity = %wd\n\n", rank, nullity);
                    flint_abort();
                }

                gr_mat_init(AX, r, nullity, ctx);
                status2 |= gr_mat_mul(AX, A, X, ctx);

                if (status2 == GR_SUCCESS && gr_mat_is_zero(AX, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (2):\n");
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); gr_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("AX: "); gr_mat_print(AX, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, nullity = %wd\n\n", rank, nullity);
                    flint_abort();
                }

                status2 |= gr_mat_rank(&nullity, X, ctx);

                if (status2 == GR_SUCCESS && nullity != gr_mat_ncols(X, ctx))
                {
                    flint_printf("FAIL (3):\n");
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); gr_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank, nullity);
                    flint_abort();
                }

                gr_mat_clear(AX, ctx);
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(X, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
