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
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_GR_FUNCTION_START(gr_mat_exp, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong n;
        gr_ctx_t ctx;
        gr_mat_t A, B, C;
        int status = GR_SUCCESS;
        int which;

        while (1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_field(ctx) == T_TRUE)
                break;
            gr_ctx_clear(ctx);
        }

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 3);
        else
            n = n_randint(state, 6);

        which = n_randint(state, 4);

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C, n, n, ctx);

        if (n_randint(state, 2))
        {
            GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        }
        else
        {
            fmpq_mat_t Q;
            fmpq_mat_init(Q, n, n);
            fmpq_mat_randtest(Q, state, 5);
            status = gr_mat_set_fmpq_mat(A, Q, ctx);
            fmpq_mat_clear(Q);
        }

        if (which & 1)
            status = gr_mat_log(B, A, ctx);
        else
            status = gr_mat_log_jordan(B, A, ctx);

        if (status == GR_SUCCESS)
        {
            if (which & 3)
                status = gr_mat_exp_jordan(C, B, ctx);
            else
                status = gr_mat_exp(C, B, ctx);
        }

        if (status == GR_SUCCESS && gr_mat_equal(A, C, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "), gr_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "), gr_mat_print(C, ctx); flint_printf("\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
