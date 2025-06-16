/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_permanent, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status = GR_SUCCESS;
        ulong m;
        slong n;
        gr_ctx_t ctx;
        gr_mat_t A;
        gr_ptr a, b;
        int algorithm;

        switch (n_randint(state, 3) % 3)
        {
            case 0:
                gr_ctx_init_fmpz(ctx);
                n = n_randint(state, 7);
                break;
            case 1:
                gr_ctx_init_fmpq(ctx);
                n = n_randint(state, 3);
                break;
            default:
                m = n_randtest_not_zero(state);
                gr_ctx_init_nmod(ctx, m);
                n = n_randint(state, 7);
                break;
        }

        gr_mat_init(A, n, n, ctx);

        a = gr_heap_init(ctx);
        b = gr_heap_init(ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

        status |= gr_mat_permanent_cofactor(a, A, ctx);

        for (algorithm = 0; algorithm < 4; algorithm++)
        {
            status = GR_SUCCESS;

            switch (algorithm)
            {
                case 0:
                    status |= gr_mat_permanent_ryser(b, A, ctx);
                    break;
                case 1:
                    status |= gr_mat_permanent_glynn(b, A, ctx);
                    break;
                case 2:
                    flint_set_num_threads(1 + n_randint(state, 5));
                    status |= gr_mat_permanent_glynn_threaded(b, A, ctx);
                    break;
                default:
                    status |= gr_mat_permanent(b, A, ctx);
                    break;
            }

            if (status == GR_SUCCESS && gr_equal(a, b, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n");
                flint_printf("algorithm = %d\n", algorithm);
                gr_ctx_println(ctx);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("a = "), gr_println(a, ctx);
                flint_printf("b = "), gr_println(b, ctx);
                flint_abort();
            }
        }

        gr_mat_clear(A, ctx);

        gr_heap_clear(a, ctx);
        gr_heap_clear(b, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
