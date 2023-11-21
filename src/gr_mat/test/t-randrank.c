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

TEST_FUNCTION_START(gr_mat_randrank, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        slong rank1, rank2, r, c;
        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        if (gr_ctx_is_integral_domain(ctx) == T_TRUE)
        {
            r = n_randint(state, 6);
            c = n_randint(state, 6);

            gr_mat_init(A, r, c, ctx);

            rank1 = n_randint(state, 6);

            status |= gr_mat_randrank(A, state, rank1, ctx);
            status |= gr_mat_rank(&rank2, A, ctx);

            if (status == GR_SUCCESS)
            {
                if (rank1 != rank2)
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank1, rank2);
                    flint_abort();
                }
            }

            gr_mat_clear(A, ctx);
        }

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
