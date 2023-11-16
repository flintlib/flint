/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_mat.h"

TEST_FUNCTION_START(ca_mat_rank, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B;
        slong rank1, rank2, r, c;
        int success;

        ca_ctx_init(ctx);

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        ca_mat_init(A, r, c, ctx);
        ca_mat_init(B, r, c, ctx);

        ca_mat_randtest(A, state, 1, 5, ctx);
        ca_mat_set(B, A, ctx);
        ca_mat_randops(B, state, 1 + n_randint(state, 20), ctx);

        success = ca_mat_rank(&rank1, A, ctx);

        if (success)
        {
            success = ca_mat_rank(&rank2, B, ctx);

            if (success)
            {
                if (rank1 != rank2)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); ca_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank1, rank2);
                    flint_abort();
                }
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
