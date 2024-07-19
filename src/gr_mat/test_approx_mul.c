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

void gr_mat_test_approx_mul_max_norm(gr_method_mat_binary_op mul_impl, gr_srcptr rel_tol, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_mat_t A, B, C, D, ERR;
        gr_ptr err, amag, bmag, tol;
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
        gr_mat_init(ERR, a, c, ctx);
        err = gr_heap_init(ctx);
        amag = gr_heap_init(ctx);
        bmag = gr_heap_init(ctx);
        tol = gr_heap_init(ctx);

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

        status |= gr_mat_sub(ERR, C, D, ctx);

        status |= gr_mat_norm_max(err, ERR, ctx);
        status |= gr_mat_norm_max(amag, A, ctx);
        status |= gr_mat_norm_max(bmag, B, ctx);
        status |= gr_mul(tol, amag, bmag, ctx);
        status |= gr_mul(tol, tol, rel_tol, ctx);

        if (status == GR_SUCCESS && gr_le(err, tol, ctx) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            gr_ctx_println(ctx);
            flint_printf("A:\n"); gr_mat_print(A, ctx); flint_printf("\n\n");
            flint_printf("B:\n"); gr_mat_print(B, ctx); flint_printf("\n\n");
            flint_printf("C:\n"); gr_mat_print(C, ctx); flint_printf("\n\n");
            flint_printf("D:\n"); gr_mat_print(D, ctx); flint_printf("\n\n");
            flint_printf("ERR:\n"); gr_mat_print(ERR, ctx); flint_printf("\n\n");
            flint_printf("err:\n"); gr_println(err, ctx); flint_printf("\n\n");
            flint_printf("tol:\n"); gr_println(tol, ctx); flint_printf("\n\n");
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);
        gr_heap_clear(err, ctx);
        gr_heap_clear(amag, ctx);
        gr_heap_clear(bmag, ctx);
        gr_heap_clear(tol, ctx);

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}

void gr_mat_test_approx_mul_pos_entrywise_accurate(gr_method_mat_binary_op mul_impl, gr_srcptr rel_tol, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_mat_t A, B, C, D, ERR, TOL;
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
        gr_mat_init(ERR, a, c, ctx);
        gr_mat_init(TOL, a, c, ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_randtest(B, state, ctx);
        status |= gr_mat_entrywise_unary_op(A, (gr_method_unary_op) gr_abs, A, ctx);
        status |= gr_mat_entrywise_unary_op(B, (gr_method_unary_op) gr_abs, B, ctx);

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

        /* |C-D| <= |D| tol */
        status |= gr_mat_sub(ERR, C, D, ctx);
        status |= gr_mat_entrywise_unary_op(ERR, (gr_method_unary_op) gr_abs, ERR, ctx);
        status |= gr_mat_entrywise_unary_op(TOL, (gr_method_unary_op) gr_abs, D, ctx);
        status |= gr_mat_mul_scalar(TOL, TOL, rel_tol, ctx);

        if (status == GR_SUCCESS && gr_mat_entrywise_binary_predicate_all((gr_method_binary_predicate) gr_le, ERR, TOL, ctx) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            gr_ctx_println(ctx);
            flint_printf("A:\n"); gr_mat_print(A, ctx); flint_printf("\n\n");
            flint_printf("B:\n"); gr_mat_print(B, ctx); flint_printf("\n\n");
            flint_printf("C:\n"); gr_mat_print(C, ctx); flint_printf("\n\n");
            flint_printf("D:\n"); gr_mat_print(D, ctx); flint_printf("\n\n");
            flint_printf("ERR:\n"); gr_mat_print(ERR, ctx); flint_printf("\n\n");
            flint_printf("TOL:\n"); gr_mat_print(TOL, ctx); flint_printf("\n\n");
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);
        gr_mat_clear(ERR, ctx);
        gr_mat_clear(TOL, ctx);

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}
