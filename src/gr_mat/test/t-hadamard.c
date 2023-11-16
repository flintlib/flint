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

TEST_GR_FUNCTION_START(gr_mat_hadamard, state, count_success, count_unable, count_domain)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B, C;
        slong n;
        int status = GR_SUCCESS;

        if (n_randint(state, 20) == 0)
        {
            gr_ctx_init_fmpz(ctx);
            n = 4 * n_randint(state, 30);
        }
        else
        {
            gr_ctx_init_random(ctx, state);
            n = n_randint(state, 20);
        }

        gr_mat_init(A, n, n, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        status = gr_mat_hadamard(A, ctx);

        if (status == GR_SUCCESS)
        {
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            status |= gr_mat_transpose(B, A, ctx);
            status |= gr_mat_mul(C, A, B, ctx);
            status |= gr_mat_set_ui(B, n, ctx);

            if (status == GR_SUCCESS && gr_mat_equal(C, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL:\n");
                gr_ctx_println(ctx);
                flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B: "); gr_mat_print(B, ctx); flint_printf("\n");
                flint_printf("C: "); gr_mat_print(C, ctx); flint_printf("\n");
                flint_abort();
            }

            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
