/*
    Copyright (C) 2022, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you grn redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

static void _gr_mat_test_nonsingular_solve_tri(gr_method_mat_binary_op_with_flag solve_impl, int upper, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t my_ctx;
        gr_ctx_ptr ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        {
            gr_mat_t A, X, B, Y;
            slong rows, cols, i, j;
            int unit;
            int status = GR_SUCCESS;
            slong sz;

            rows = n_randint(state, maxn + 1);
            cols = n_randint(state, maxn + 1);
            unit = n_randint(state, 2);

            sz = ((gr_ctx_struct *) ctx)->sizeof_elem;

            gr_mat_init(A, rows, rows, ctx);
            gr_mat_init(B, rows, cols, ctx);
            gr_mat_init(X, rows, cols, ctx);
            gr_mat_init(Y, rows, cols, ctx);

            status |= gr_mat_randtest(A, state, ctx);
            status |= gr_mat_randtest(X, state, ctx);
            status |= gr_mat_randtest(Y, state, ctx);

            for (i = 0; i < rows; i++)
            {
                if (unit)
                    status |= gr_one(GR_MAT_ENTRY(A, i, i, sz), ctx);
                else
                    status |= gr_set_ui(GR_MAT_ENTRY(A, i, i, sz), 1 + n_randint(state, 100), ctx);

                if (upper)
                    for (j = 0; j < i; j++)
                        status |= gr_zero(GR_MAT_ENTRY(A, i, j, sz), ctx);
                else
                    for (j = i + 1; j < rows; j++)
                        status |= gr_zero(GR_MAT_ENTRY(A, i, j, sz), ctx);
            }

            status |= gr_mat_mul(B, A, X, ctx);

            if (unit)  /* check that diagonal entries are ignored */
            {
                for (i = 0; i < rows; i++)
                    status |= gr_set_ui(GR_MAT_ENTRY(A, i, i, sz), 1 + n_randint(state, 100), ctx);
            }

            /* Check Y = A^(-1) * (A * X) = X */
            if (n_randint(state, 2))
            {
                if (upper)
                    status |= gr_mat_nonsingular_solve_triu(Y, A, B, unit, ctx);
                else
                    status |= gr_mat_nonsingular_solve_tril(Y, A, B, unit, ctx);
            }
            else
            {
                status |= gr_mat_set(Y, B, ctx);
                if (upper)
                    status |= gr_mat_nonsingular_solve_triu(Y, A, Y, unit, ctx);
                else
                    status |= gr_mat_nonsingular_solve_tril(Y, A, Y, unit, ctx);
            }

            if (status == GR_SUCCESS && gr_mat_equal(Y, X, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n");
                gr_ctx_println(ctx);
                flint_printf("A = \n"); gr_mat_print(A, ctx); flint_printf("\n\n");
                flint_printf("B = \n"); gr_mat_print(B, ctx); flint_printf("\n\n");
                flint_printf("X = \n"); gr_mat_print(X, ctx); flint_printf("\n\n");
                flint_printf("Y = \n"); gr_mat_print(Y, ctx); flint_printf("\n\n");
                flint_abort();
            }

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(X, ctx);
            gr_mat_clear(Y, ctx);
        }

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}

void gr_mat_test_nonsingular_solve_tril(gr_method_mat_binary_op_with_flag solve_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    _gr_mat_test_nonsingular_solve_tri(solve_impl, 0, state, iters, maxn, ctx);
}

void gr_mat_test_nonsingular_solve_triu(gr_method_mat_binary_op_with_flag solve_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    _gr_mat_test_nonsingular_solve_tri(solve_impl, 1, state, iters, maxn, ctx);
}
