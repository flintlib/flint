/*
    Copyright (C) 2022, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you grn redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"

void gr_mat_test_mul(gr_method_mat_binary_op mul_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_mat_t A, B, C, D;
        slong a, b, c;
        int status = GR_SUCCESS;
        gr_ctx_t my_ctx;
        gr_ctx_ptr ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        if (n_randint(state, 4) == 0)
        {
            a = b = c = n_randint(state, maxn);
        }
        else
        {
            a = n_randint(state, maxn);
            b = n_randint(state, maxn);
            c = n_randint(state, maxn);
        }

        gr_mat_init(A, a, b, ctx);
        gr_mat_init(B, b, c, ctx);
        gr_mat_init(C, a, c, ctx);
        gr_mat_init(D, a, c, ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_randtest(B, state, ctx);
        status |= gr_mat_randtest(C, state, ctx);
        status |= gr_mat_randtest(D, state, ctx);

        if (b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(C, A, ctx);
            status |= mul_impl(C, C, B, ctx);
        }
        else if (a == b && n_randint(state, 2))
        {
            status |= gr_mat_set(C, B, ctx);
            status |= mul_impl(C, A, C, ctx);
        }
        else if (a == b && b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(B, A, ctx);
            status |= mul_impl(C, A, A, ctx);
        }
        else if (a == b && b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(B, A, ctx);
            status |= gr_mat_set(C, A, ctx);
            status |= mul_impl(C, C, C, ctx);
        }
        else
        {
            status |= mul_impl(C, A, B, ctx);
        }

        status |= gr_mat_mul_classical(D, A, B, ctx);

        if (status == GR_SUCCESS && gr_mat_equal(C, D, ctx) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            gr_ctx_println(ctx);
            flint_printf("A:\n"); gr_mat_print(A, ctx); flint_printf("\n\n");
            flint_printf("B:\n"); gr_mat_print(B, ctx); flint_printf("\n\n");
            flint_printf("C:\n"); gr_mat_print(C, ctx); flint_printf("\n\n");
            flint_printf("D:\n"); gr_mat_print(D, ctx); flint_printf("\n\n");
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}
