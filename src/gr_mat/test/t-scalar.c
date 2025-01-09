/*
    Copyright (C) 2022, 2025 Fredrik Johansson

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

TEST_FUNCTION_START(gr_mat_scalar, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx, ctx_other;
        gr_mat_t A, B1, B2, B3, Cmn, Cmm, Cnn;
        gr_ptr c, d, c_other;
        slong m, n;
        int status = GR_SUCCESS;
        int testcase;
        int have_other;

        m = n_randint(state, 4);
        n = n_randint(state, 4);

        /* Hack: avoid because slow */
        gr_ctx_init_random(ctx, state);
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        gr_ctx_init_random(ctx_other, state);
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx_other);
            gr_ctx_init_random(ctx_other, state);
        }

        gr_mat_init(A, m, n, ctx);
        gr_mat_init(Cmn, m, n, ctx);
        gr_mat_init(Cmm, m, m, ctx);
        gr_mat_init(Cnn, n, n, ctx);
        gr_mat_init(B1, m, n, ctx);
        gr_mat_init(B2, m, n, ctx);
        gr_mat_init(B3, m, n, ctx);
        c = gr_heap_init(ctx);
        d = gr_heap_init(ctx);
        c_other = gr_heap_init(ctx_other);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        GR_MUST_SUCCEED(gr_randtest(c, state, ctx));

        have_other = (gr_set_other(c_other, c, ctx, ctx_other) == GR_SUCCESS);

        status |= gr_mat_set_scalar(Cmn, c, ctx);
        status |= gr_mat_set_scalar(Cmm, c, ctx);
        status |= gr_mat_set_scalar(Cnn, c, ctx);

        for (testcase = 0; testcase <= 6; testcase++)
        {
            status = GR_SUCCESS;

            gr_mat_struct * A_or_B2_alias;
            gr_mat_struct * A_or_B3_alias;

            if (n_randint(state, 2))
            {
                A_or_B2_alias = A;
                A_or_B3_alias = A;
            }
            else
            {
                status |= gr_mat_set(B2, A, ctx);
                status |= gr_mat_set(B3, A, ctx);
                A_or_B2_alias = B2;
                A_or_B3_alias = B3;
            }

            if (testcase == 0)
            {
                /* A + c == A + C */
                status |= gr_mat_add(B1, A, Cmn, ctx);
                status |= gr_mat_add_scalar(B2, A_or_B2_alias, c, ctx);
                if (have_other)
                    status |= gr_mat_add_scalar_other(B3, A_or_B3_alias, c_other, ctx_other, ctx);
            }
            else if (testcase == 1)
            {
                /* c + A == C + A */
                status |= gr_mat_add(B1, Cmn, A, ctx);
                status |= gr_mat_scalar_add(B2, c, A_or_B2_alias, ctx);
                if (have_other)
                    status |= gr_mat_scalar_other_add(B3, c_other, ctx_other, A_or_B3_alias, ctx);
            }
            else if (testcase == 2)
            {
                /* A - c == A - C */
                status |= gr_mat_sub(B1, A, Cmn, ctx);
                status |= gr_mat_sub_scalar(B2, A_or_B2_alias, c, ctx);
                if (have_other)
                    status |= gr_mat_sub_scalar_other(B3, A_or_B3_alias, c_other, ctx_other, ctx);
            }
            else if (testcase == 3)
            {
                /* c - A == C - A */
                status |= gr_mat_sub(B1, Cmn, A, ctx);
                status |= gr_mat_scalar_sub(B2, c, A_or_B2_alias, ctx);
                if (have_other)
                    status |= gr_mat_scalar_other_sub(B3, c_other, ctx_other, A_or_B3_alias, ctx);
            }
            else if (testcase == 4)
            {
                /* A * c == A * C */
                status |= gr_mat_mul(B1, A, Cnn, ctx);
                status |= gr_mat_mul_scalar(B2, A_or_B2_alias, c, ctx);
                if (have_other)
                    status |= gr_mat_mul_scalar_other(B3, A_or_B3_alias, c_other, ctx_other, ctx);
            }
            else if (testcase == 5)
            {
                /* A * c == A * C */
                status |= gr_mat_mul(B1, Cmm, A, ctx);
                status |= gr_mat_scalar_mul(B2, c, A_or_B2_alias, ctx);
                if (have_other)
                    status |= gr_mat_scalar_other_mul(B3, c_other, ctx_other, A_or_B3_alias, ctx);
            }
            else if (testcase == 6)
            {
                /* A / c == A * c^(-1) */
                status |= gr_inv(d, c, ctx);
                status |= gr_mat_mul_scalar(B1, A, d, ctx);
                status |= gr_mat_div_scalar(B2, A_or_B2_alias, c, ctx);
                if (have_other)
                    status |= gr_mat_div_scalar_other(B3, A_or_B3_alias, c_other, ctx_other, ctx);
            }

            if (status == GR_SUCCESS && gr_mat_equal(B1, B2, ctx) == T_FALSE)
            {
                flint_printf("FAIL (scalar %d)\n", testcase);
                gr_ctx_println(ctx);
                flint_printf("A = \n"); gr_mat_print(A, ctx); flint_printf("\n\n");
                flint_printf("c = \n"); gr_print(c, ctx); flint_printf("\n\n");
                flint_printf("Cmn = \n"); gr_mat_print(Cmn, ctx); flint_printf("\n\n");
                flint_printf("Cmm = \n"); gr_mat_print(Cmm, ctx); flint_printf("\n\n");
                flint_printf("Cnn = \n"); gr_mat_print(Cnn, ctx); flint_printf("\n\n");
                flint_printf("B1 = \n"); gr_mat_print(B1, ctx); flint_printf("\n\n");
                flint_printf("B2 = \n"); gr_mat_print(B2, ctx); flint_printf("\n\n");
                flint_abort();
            }

            if (status == GR_SUCCESS && have_other && gr_mat_equal(B1, B3, ctx) == T_FALSE)
            {
                flint_printf("FAIL (scalar %d)\n", testcase);
                gr_ctx_println(ctx);
                flint_printf("A = \n"); gr_mat_print(A, ctx); flint_printf("\n\n");
                flint_printf("c = \n"); gr_print(c, ctx); flint_printf("\n\n");
                flint_printf("Cmn = \n"); gr_mat_print(Cmn, ctx); flint_printf("\n\n");
                flint_printf("Cmm = \n"); gr_mat_print(Cmm, ctx); flint_printf("\n\n");
                flint_printf("Cnn = \n"); gr_mat_print(Cnn, ctx); flint_printf("\n\n");
                flint_printf("B1 = \n"); gr_mat_print(B1, ctx); flint_printf("\n\n");
                flint_printf("B3 = \n"); gr_mat_print(B3, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(Cmn, ctx);
        gr_mat_clear(Cmm, ctx);
        gr_mat_clear(Cnn, ctx);
        gr_mat_clear(B1, ctx);
        gr_mat_clear(B2, ctx);
        gr_mat_clear(B3, ctx);
        gr_heap_clear(c, ctx);
        gr_heap_clear(d, ctx);
        gr_heap_clear(c_other, ctx_other);

        gr_ctx_clear(ctx);
        gr_ctx_clear(ctx_other);
    }

    TEST_FUNCTION_END(state);
}
