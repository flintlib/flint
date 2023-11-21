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

TEST_FUNCTION_START(gr_mat_properties, state)
{
    slong iter;

    for (iter = 0; iter < 10000; iter++)
    {
        int status = GR_SUCCESS;
        slong i, j, m, n;
        gr_ctx_t ctx;
        gr_mat_t A, B;
        truth_t eq1, eq2;

        gr_ctx_init_fmpz(ctx);

        m = n_randint(state, 5);
        n = n_randint(state, 5);
        gr_mat_init(A, m, n, ctx);
        gr_mat_init(B, m, n, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

        eq1 = gr_mat_is_lower_triangular(A, ctx);
        status = gr_mat_set(B, A, ctx);
        for (i = 0; i < m; i++)
            for (j = i + 1; j < n; j++)
                status |= gr_zero(gr_mat_entry_ptr(B, i, j, ctx), ctx);
        eq2 = gr_mat_equal(A, B, ctx);

        if (status == GR_SUCCESS && ((eq1 == T_TRUE && eq2 == T_FALSE) || (eq1 == T_FALSE && eq2 == T_TRUE)))
        {
            flint_printf("FAIL (is_lower_triangular)\n");
            printf("A = "); gr_mat_print(A, ctx); printf("\n");
            printf("B = "); gr_mat_print(B, ctx); printf("\n");
            truth_println(eq1);
            truth_println(eq2);
            flint_abort();
        }

        eq1 = gr_mat_is_upper_triangular(A, ctx);
        status = gr_mat_set(B, A, ctx);
        for (i = 0; i < m; i++)
            for (j = 0; j < FLINT_MIN(i, n); j++)
                status |= gr_zero(gr_mat_entry_ptr(B, i, j, ctx), ctx);
        eq2 = gr_mat_equal(A, B, ctx);

        if (status == GR_SUCCESS && ((eq1 == T_TRUE && eq2 == T_FALSE) || (eq1 == T_FALSE && eq2 == T_TRUE)))
        {
            flint_printf("FAIL (is_upper_triangular)\n");
            printf("A = "); gr_mat_print(A, ctx); printf("\n");
            printf("B = "); gr_mat_print(B, ctx); printf("\n");
            truth_println(eq1);
            truth_println(eq2);
            flint_abort();
        }

        eq1 = gr_mat_is_hessenberg(A, ctx);
        status = gr_mat_set(B, A, ctx);
        for (i = 0; i < m; i++)
            for (j = 0; j < FLINT_MIN(i - 1, n); j++)
                status |= gr_zero(gr_mat_entry_ptr(B, i, j, ctx), ctx);
        eq2 = gr_mat_equal(A, B, ctx);

        if (status == GR_SUCCESS && ((eq1 == T_TRUE && eq2 == T_FALSE) || (eq1 == T_FALSE && eq2 == T_TRUE)))
        {
            flint_printf("FAIL (is_hessenberg)\n");
            printf("A = "); gr_mat_print(A, ctx); printf("\n");
            printf("B = "); gr_mat_print(B, ctx); printf("\n");
            truth_println(eq1);
            truth_println(eq2);
            flint_abort();
        }

        eq1 = gr_mat_is_diagonal(A, ctx);
        status = gr_mat_set(B, A, ctx);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                if (i != j)
                    status |= gr_zero(gr_mat_entry_ptr(B, i, j, ctx), ctx);
        eq2 = gr_mat_equal(A, B, ctx);

        if (status == GR_SUCCESS && ((eq1 == T_TRUE && eq2 == T_FALSE) || (eq1 == T_FALSE && eq2 == T_TRUE)))
        {
            flint_printf("FAIL (is_diagonal)\n");
            printf("A = "); gr_mat_print(A, ctx); printf("\n");
            printf("B = "); gr_mat_print(B, ctx); printf("\n");
            truth_println(eq1);
            truth_println(eq2);
            flint_abort();
        }

        eq1 = gr_mat_is_zero(A, ctx);
        status = gr_mat_set(B, A, ctx);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                status |= gr_zero(gr_mat_entry_ptr(B, i, j, ctx), ctx);
        eq2 = gr_mat_equal(A, B, ctx);

        if (status == GR_SUCCESS && ((eq1 == T_TRUE && eq2 == T_FALSE) || (eq1 == T_FALSE && eq2 == T_TRUE)))
        {
            flint_printf("FAIL (is_zero)\n");
            printf("A = "); gr_mat_print(A, ctx); printf("\n");
            printf("B = "); gr_mat_print(B, ctx); printf("\n");
            truth_println(eq1);
            truth_println(eq2);
            flint_abort();
        }

        eq1 = gr_mat_is_one(A, ctx);
        status = gr_mat_set(B, A, ctx);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                if (i == j)
                    status |= gr_one(gr_mat_entry_ptr(B, i, j, ctx), ctx);
                else
                    status |= gr_zero(gr_mat_entry_ptr(B, i, j, ctx), ctx);
        eq2 = gr_mat_equal(A, B, ctx);

        if (status == GR_SUCCESS && ((eq1 == T_TRUE && eq2 == T_FALSE) || (eq1 == T_FALSE && eq2 == T_TRUE)))
        {
            flint_printf("FAIL (is_one)\n");
            printf("A = "); gr_mat_print(A, ctx); printf("\n");
            printf("B = "); gr_mat_print(B, ctx); printf("\n");
            truth_println(eq1);
            truth_println(eq2);
            flint_abort();
        }

        eq1 = gr_mat_is_neg_one(A, ctx);
        status = gr_mat_set(B, A, ctx);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                if (i == j)
                    status |= gr_neg_one(gr_mat_entry_ptr(B, i, j, ctx), ctx);
                else
                    status |= gr_zero(gr_mat_entry_ptr(B, i, j, ctx), ctx);
        eq2 = gr_mat_equal(A, B, ctx);

        if (status == GR_SUCCESS && ((eq1 == T_TRUE && eq2 == T_FALSE) || (eq1 == T_FALSE && eq2 == T_TRUE)))
        {
            flint_printf("FAIL (is_neg_one)\n");
            printf("A = "); gr_mat_print(A, ctx); printf("\n");
            printf("B = "); gr_mat_print(B, ctx); printf("\n");
            truth_println(eq1);
            truth_println(eq2);
            flint_abort();
        }

        eq1 = gr_mat_is_scalar(A, ctx);
        status = gr_mat_set(B, A, ctx);
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
                if (i != j)
                    status |= gr_zero(gr_mat_entry_ptr(B, i, j, ctx), ctx);
            if (i > 0 && i < n)
                status |= gr_set(gr_mat_entry_ptr(B, i, i, ctx), gr_mat_entry_ptr(B, 0, 0, ctx), ctx);
        }
        eq2 = gr_mat_equal(A, B, ctx);

        if (status == GR_SUCCESS && ((eq1 == T_TRUE && eq2 == T_FALSE) || (eq1 == T_FALSE && eq2 == T_TRUE)))
        {
            flint_printf("FAIL (is_scalar)\n");
            printf("A = "); gr_mat_print(A, ctx); printf("\n");
            printf("B = "); gr_mat_print(B, ctx); printf("\n");
            truth_println(eq1);
            truth_println(eq2);
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
