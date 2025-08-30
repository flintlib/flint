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

TEST_GR_FUNCTION_START(gr_mat_qr, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        slong n, m;
        gr_ctx_t ctx;
        gr_mat_t A, R, Q, QR;
        int alias = n_randint(state, 3);

        gr_ctx_init_random(ctx, state);

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 3);
        else
            n = n_randint(state, 8);

        m = n_randint(state, 1 + n);

        if (m != n && alias == 2)
            alias = 1;

        gr_mat_init(A, n, m, ctx);
        gr_mat_init(R, m, m, ctx);
        gr_mat_init(Q, n, m, ctx);
        gr_mat_init(QR, n, m, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(R, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(Q, state, ctx));

        if (alias == 2)
        {
            status |= gr_mat_set(R, A, ctx);
            status |= gr_mat_qr(Q, R, R, ctx);
        }
        else if (alias == 1)
        {
            status |= gr_mat_set(Q, A, ctx);
            status |= gr_mat_qr(Q, R, Q, ctx);
        }
        else
        {
            status |= gr_mat_qr(Q, R, A, ctx);
        }

        if (status == GR_SUCCESS)
        {
            if (gr_mat_is_upper_triangular(R, ctx) == T_FALSE)
            {
                flint_printf("FAIL: R not lower triangular\n");
                gr_ctx_println(ctx);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("Q = "), gr_mat_print(Q, ctx); flint_printf("\n");
                flint_printf("R = "), gr_mat_print(R, ctx); flint_printf("\n");
                flint_abort();
            }

            if (gr_mat_is_col_orthogonal(Q, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Q not orthogonal\n");
                gr_ctx_println(ctx);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("Q = "), gr_mat_print(Q, ctx); flint_printf("\n");
                flint_printf("R = "), gr_mat_print(R, ctx); flint_printf("\n");
                flint_abort();
            }

            status |= gr_mat_mul(QR, Q, R, ctx);

            if (status == GR_SUCCESS && gr_mat_equal(QR, A, ctx) == T_FALSE)
            {
                flint_printf("FAIL: QR != A\n");
                gr_ctx_println(ctx);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("Q = "), gr_mat_print(Q, ctx); flint_printf("\n");
                flint_printf("R = "), gr_mat_print(R, ctx); flint_printf("\n");
                flint_printf("Q*R = "), gr_mat_print(QR, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(R, ctx);
        gr_mat_clear(Q, ctx);
        gr_mat_clear(QR, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
