/*
    Copyright (C) 2022 Fredrik Johansson

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

TEST_GR_FUNCTION_START(gr_mat_lq, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        slong n, m;
        gr_ctx_t ctx;
        gr_mat_t A, L, Q, LQ;
        int alias = n_randint(state, 3);
        int method = n_randint(state, 3);
        gr_method_mat_binary_unary_op func;

        gr_ctx_init_random(ctx, state);

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 3);
        else
            n = n_randint(state, 8);

        m = n_randint(state, 1 + n);

        if (m != n && alias == 2)
            alias = 1;

        gr_mat_init(A, m, n, ctx);
        gr_mat_init(L, m, m, ctx);
        gr_mat_init(Q, m, n, ctx);
        gr_mat_init(LQ, m, n, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(L, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(Q, state, ctx));

        if (method == 0)
            func = (gr_method_mat_binary_unary_op) gr_mat_lq;
        else if (method == 1)
            func = (gr_method_mat_binary_unary_op) gr_mat_lq_gso;
        else
            func = (gr_method_mat_binary_unary_op) gr_mat_lq_recursive;

        if (alias == 2)
        {
            status |= gr_mat_set(L, A, ctx);
            status |= func(L, Q, L, ctx);
        }
        else if (alias == 1)
        {
            status |= gr_mat_set(Q, A, ctx);
            status |= func(L, Q, Q, ctx);
        }
        else
        {
            status |= func(L, Q, A, ctx);
        }

        if (status == GR_SUCCESS)
        {
            if (gr_mat_is_lower_triangular(L, ctx) == T_FALSE)
            {
                flint_printf("FAIL: L not lower triangular\n");
                gr_ctx_println(ctx);
                flint_printf("alias = %wd, method = %wd\n", alias, method);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("L = "), gr_mat_print(L, ctx); flint_printf("\n");
                flint_printf("Q = "), gr_mat_print(Q, ctx); flint_printf("\n");
                flint_abort();
            }

            if (gr_mat_is_row_orthogonal(Q, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Q not orthogonal\n");
                gr_ctx_println(ctx);
                flint_printf("alias = %wd, method = %wd\n", alias, method);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("L = "), gr_mat_print(L, ctx); flint_printf("\n");
                flint_printf("Q = "), gr_mat_print(Q, ctx); flint_printf("\n");
                flint_abort();
            }

            status |= gr_mat_mul(LQ, L, Q, ctx);

            if (status == GR_SUCCESS && gr_mat_equal(LQ, A, ctx) == T_FALSE)
            {
                flint_printf("FAIL: LQ != A\n");
                gr_ctx_println(ctx);
                flint_printf("alias = %wd, method = %wd\n", alias, method);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("L = "), gr_mat_print(L, ctx); flint_printf("\n");
                flint_printf("Q = "), gr_mat_print(Q, ctx); flint_printf("\n");
                flint_printf("L*Q = "), gr_mat_print(LQ, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(L, ctx);
        gr_mat_clear(Q, ctx);
        gr_mat_clear(LQ, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
