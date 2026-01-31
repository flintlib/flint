/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_mat.h"

static truth_t
gr_mat_is_orthogonal_naive(const gr_mat_t A, gr_ctx_t ctx)
{
    gr_mat_t T, U;
    slong n = A->r;
    truth_t res, res2;
    int status = GR_SUCCESS;

    gr_mat_init(T, n, n, ctx);
    gr_mat_init(U, n, n, ctx);

    status |= gr_mat_transpose(T, A, ctx);
    status |= gr_mat_mul(U, A, T, ctx);
    res = gr_mat_is_one(U, ctx);

    status |= gr_mat_mul(U, T, A, ctx);
    res2 = gr_mat_is_one(U, ctx);
    res = truth_and(res, res2);

    gr_mat_clear(T, ctx);
    gr_mat_clear(U, ctx);

    if (status != GR_SUCCESS)
        res = T_UNKNOWN;

    return res;
}

static truth_t
gr_mat_is_orthogonal2_naive(const gr_mat_t A, int cols, int unit, gr_ctx_t ctx)
{
    gr_mat_t AT, P;
    slong r = A->r;
    slong c = A->c;
    truth_t res;
    int status = GR_SUCCESS;

    gr_mat_init(AT, c, r, ctx);
    status |= gr_mat_transpose(AT, A, ctx);

    if (cols)
    {
        gr_mat_init(P, c, c, ctx);
        status |= gr_mat_mul(P, AT, A, ctx);
    }
    else
    {
        gr_mat_init(P, r, r, ctx);
        status |= gr_mat_mul(P, A, AT, ctx);
    }

    if (unit)
        res = gr_mat_is_one(P, ctx);
    else
        res = gr_mat_is_diagonal(P, ctx);

    gr_mat_clear(AT, ctx);
    gr_mat_clear(P, ctx);

    if (status != GR_SUCCESS)
        res = T_UNKNOWN;

    return res;
}

TEST_FUNCTION_START(gr_mat_is_orthogonal, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong r, c;
        gr_ctx_t ctx;
        gr_mat_t A;
        truth_t r1, r2;

        gr_ctx_init_nmod(ctx, n_randtest_not_zero(state));

        r = n_randint(state, 10);
        c = n_randint(state, 10);

        gr_mat_init(A, r, c, ctx);
        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

        r1 = gr_mat_is_row_orthogonal(A, ctx);
        r2 = gr_mat_is_orthogonal2_naive(A, 0, 0, ctx);

        if ((r1 == T_FALSE && r2 == T_TRUE) || (r1 == T_TRUE && r2 == T_FALSE))
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("is_row_orthogonal\n");
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("r1 = %{truth}\n", r1);
            flint_printf("r2 = %{truth}\n", r2);
            flint_abort();
        }

        r1 = gr_mat_is_row_orthonormal(A, ctx);
        r2 = gr_mat_is_orthogonal2_naive(A, 0, 1, ctx);

        if ((r1 == T_FALSE && r2 == T_TRUE) || (r1 == T_TRUE && r2 == T_FALSE))
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("is_row_orthonormal\n");
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("r1 = %{truth}\n", r1);
            flint_printf("r2 = %{truth}\n", r2);
            flint_abort();
        }

        r1 = gr_mat_is_col_orthogonal(A, ctx);
        r2 = gr_mat_is_orthogonal2_naive(A, 1, 0, ctx);

        if ((r1 == T_FALSE && r2 == T_TRUE) || (r1 == T_TRUE && r2 == T_FALSE))
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("is_col_orthonormal\n");
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("r1 = %{truth}\n", r1);
            flint_printf("r2 = %{truth}\n", r2);
            flint_abort();
        }

        r1 = gr_mat_is_col_orthonormal(A, ctx);
        r2 = gr_mat_is_orthogonal2_naive(A, 1, 1, ctx);

        if ((r1 == T_FALSE && r2 == T_TRUE) || (r1 == T_TRUE && r2 == T_FALSE))
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("is_col_orthonormal\n");
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("r1 = %{truth}\n", r1);
            flint_printf("r2 = %{truth}\n", r2);
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong n;
        gr_ctx_t ctx;
        gr_mat_t A;
        truth_t r1, r2;

        switch (n_randint(state, 4))
        {
            case 0:
                gr_ctx_init_fmpz(ctx);
                n = n_randint(state, 15);
                break;
            case 1:
                gr_ctx_init_fmpq(ctx);
                n = n_randint(state, 5);
                break;
            case 2:
                gr_ctx_init_nmod(ctx, n_randtest_not_zero(state));
                n = n_randint(state, 50);
                break;
            default:
                gr_ctx_init_random(ctx, state);
                n = n_randint(state, 4);
                break;
        }

        gr_mat_init(A, n, n, ctx);
        GR_MUST_SUCCEED(gr_mat_randtest_orthogonal(A, state, ctx));

        r1 = gr_mat_is_orthogonal(A, ctx);

        if (r1 == T_FALSE)
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("r1 = %{truth}\n", r1);
            flint_abort();
        }

        if (n != 0)
        {
            GR_IGNORE(gr_randtest(gr_mat_entry_ptr(A,
                n_randint(state, n), n_randint(state, n), ctx), state, ctx));

            r1 = gr_mat_is_orthogonal(A, ctx);
            r2 = gr_mat_is_orthogonal_naive(A, ctx);

            if ((r1 == T_FALSE && r2 == T_TRUE) || (r1 == T_TRUE && r2 == T_FALSE))
            {
                flint_printf("FAIL\n");
                gr_ctx_println(ctx);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("r1 = %{truth}\n", r1);
                flint_printf("r2 = %{truth}\n", r2);
                flint_abort();
            }
        }

        GR_IGNORE(gr_mat_randtest(A, state, ctx));

        r1 = gr_mat_is_orthogonal(A, ctx);
        r2 = gr_mat_is_orthogonal_naive(A, ctx);

        if ((r1 == T_FALSE && r2 == T_TRUE) || (r1 == T_TRUE && r2 == T_FALSE))
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("r1 = %{truth}\n", r1);
            flint_printf("r2 = %{truth}\n", r2);
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
