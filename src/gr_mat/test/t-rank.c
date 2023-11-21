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

TEST_GR_FUNCTION_START(gr_mat_rank, state, count_success, count_unable, count_domain)
{
    slong iter;

    /* Check that random row/column operations preserve rank */
    for (iter = 0; iter < 10000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B;
        slong rank1, rank2, r, c;
        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        gr_mat_init(A, r, c, ctx);
        gr_mat_init(B, r, c, ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_set(B, A, ctx);
        status |= gr_mat_randops(B, state, 1 + n_randint(state, 20), ctx);

        status |= gr_mat_rank(&rank1, A, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_mat_rank(&rank2, B, ctx);

            if (status == GR_SUCCESS)
            {
                if (rank1 != rank2)
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); gr_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank1, rank2);
                    flint_abort();
                }
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
