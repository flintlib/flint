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

    for (iter = 0; iter < iters; iter++)
    {
        gr_mat_t A, B, C, D, ERR;
        gr_ptr err, amag, bmag, tol;
        slong a, b, c;
        int status = GR_SUCCESS;
        gr_ctx_t ctx2;
        gr_ctx_ptr ctxptr;

        if (ctx == NULL)
        {
            gr_ctx_init_random(ctx2, state);
            ctxptr = ctx2;
        }
        else
            ctxptr = ctx;

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

        gr_mat_init(A, a, b, ctxptr);
        gr_mat_init(B, b, c, ctxptr);
        gr_mat_init(C, a, c, ctxptr);
        gr_mat_init(D, a, c, ctxptr);
        gr_mat_init(ERR, a, c, ctxptr);
        err = gr_heap_init(ctxptr);
        amag = gr_heap_init(ctxptr);
        bmag = gr_heap_init(ctxptr);
        tol = gr_heap_init(ctxptr);

        status |= gr_mat_randtest(A, state, ctxptr);
        status |= gr_mat_randtest(B, state, ctxptr);
        status |= gr_mat_randtest(C, state, ctxptr);
        status |= gr_mat_randtest(D, state, ctxptr);

        if (b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(C, A, ctxptr);
            status |= mul_impl(C, C, B, ctxptr);
        }
        else if (a == b && n_randint(state, 2))
        {
            status |= gr_mat_set(C, B, ctxptr);
            status |= mul_impl(C, A, C, ctxptr);
        }
        else if (a == b && b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(B, A, ctxptr);
            status |= mul_impl(C, A, A, ctxptr);
        }
        else if (a == b && b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(B, A, ctxptr);
            status |= gr_mat_set(C, A, ctxptr);
            status |= mul_impl(C, C, C, ctxptr);
        }
        else
        {
            status |= mul_impl(C, A, B, ctxptr);
        }

        status |= gr_mat_mul_classical(D, A, B, ctxptr);

        status |= gr_mat_sub(ERR, C, D, ctxptr);

        status |= gr_mat_norm_max(err, ERR, ctxptr);
        status |= gr_mat_norm_max(amag, A, ctxptr);
        status |= gr_mat_norm_max(bmag, B, ctxptr);
        status |= gr_mul(tol, amag, bmag, ctxptr);
        status |= gr_mul(tol, tol, rel_tol, ctxptr);

        if (status == GR_SUCCESS && gr_le(err, tol, ctxptr) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            gr_ctx_println(ctxptr);
            flint_printf("A:\n"); gr_mat_print(A, ctxptr); flint_printf("\n\n");
            flint_printf("B:\n"); gr_mat_print(B, ctxptr); flint_printf("\n\n");
            flint_printf("C:\n"); gr_mat_print(C, ctxptr); flint_printf("\n\n");
            flint_printf("D:\n"); gr_mat_print(D, ctxptr); flint_printf("\n\n");
            flint_printf("ERR:\n"); gr_mat_print(ERR, ctxptr); flint_printf("\n\n");
            flint_printf("err:\n"); gr_println(err, ctxptr); flint_printf("\n\n");
            flint_printf("tol:\n"); gr_println(tol, ctxptr); flint_printf("\n\n");
            flint_abort();
        }

        gr_mat_clear(A, ctxptr);
        gr_mat_clear(B, ctxptr);
        gr_mat_clear(C, ctxptr);
        gr_mat_clear(D, ctxptr);
        gr_heap_clear(err, ctxptr);
        gr_heap_clear(amag, ctxptr);
        gr_heap_clear(bmag, ctxptr);
        gr_heap_clear(tol, ctxptr);

        if (ctx == NULL)
            gr_ctx_clear(ctxptr);
    }
}

void gr_mat_test_approx_mul_pos_entrywise_accurate(gr_method_mat_binary_op mul_impl, gr_srcptr rel_tol, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        gr_mat_t A, B, C, D, ERR, TOL;
        slong a, b, c;
        int status = GR_SUCCESS;
        gr_ctx_t ctx2;
        gr_ctx_ptr ctxptr;

        if (ctx == NULL)
        {
            gr_ctx_init_random(ctx2, state);
            ctxptr = ctx2;
        }
        else
            ctxptr = ctx;

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

        gr_mat_init(A, a, b, ctxptr);
        gr_mat_init(B, b, c, ctxptr);
        gr_mat_init(C, a, c, ctxptr);
        gr_mat_init(D, a, c, ctxptr);
        gr_mat_init(ERR, a, c, ctxptr);
        gr_mat_init(TOL, a, c, ctxptr);

        status |= gr_mat_randtest(A, state, ctxptr);
        status |= gr_mat_randtest(B, state, ctxptr);
        status |= gr_mat_entrywise_unary_op(A, (gr_method_unary_op) gr_abs, A, ctxptr);
        status |= gr_mat_entrywise_unary_op(B, (gr_method_unary_op) gr_abs, B, ctxptr);

        status |= gr_mat_randtest(C, state, ctxptr);
        status |= gr_mat_randtest(D, state, ctxptr);

        if (b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(C, A, ctxptr);
            status |= mul_impl(C, C, B, ctxptr);
        }
        else if (a == b && n_randint(state, 2))
        {
            status |= gr_mat_set(C, B, ctxptr);
            status |= mul_impl(C, A, C, ctxptr);
        }
        else if (a == b && b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(B, A, ctxptr);
            status |= mul_impl(C, A, A, ctxptr);
        }
        else if (a == b && b == c && n_randint(state, 2))
        {
            status |= gr_mat_set(B, A, ctxptr);
            status |= gr_mat_set(C, A, ctxptr);
            status |= mul_impl(C, C, C, ctxptr);
        }
        else
        {
            status |= mul_impl(C, A, B, ctxptr);
        }

        status |= gr_mat_mul_classical(D, A, B, ctxptr);

        /* |C-D| <= |D| tol */
        status |= gr_mat_sub(ERR, C, D, ctxptr);
        status |= gr_mat_entrywise_unary_op(ERR, (gr_method_unary_op) gr_abs, ERR, ctxptr);
        status |= gr_mat_entrywise_unary_op(TOL, (gr_method_unary_op) gr_abs, D, ctxptr);
        status |= gr_mat_mul_scalar(TOL, TOL, rel_tol, ctxptr);

        if (status == GR_SUCCESS && gr_mat_entrywise_binary_predicate_all((gr_method_binary_predicate) gr_le, ERR, TOL, ctxptr) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            gr_ctx_println(ctxptr);
            flint_printf("A:\n"); gr_mat_print(A, ctxptr); flint_printf("\n\n");
            flint_printf("B:\n"); gr_mat_print(B, ctxptr); flint_printf("\n\n");
            flint_printf("C:\n"); gr_mat_print(C, ctxptr); flint_printf("\n\n");
            flint_printf("D:\n"); gr_mat_print(D, ctxptr); flint_printf("\n\n");
            flint_printf("ERR:\n"); gr_mat_print(ERR, ctxptr); flint_printf("\n\n");
            flint_printf("TOL:\n"); gr_mat_print(TOL, ctxptr); flint_printf("\n\n");
            flint_abort();
        }

        gr_mat_clear(A, ctxptr);
        gr_mat_clear(B, ctxptr);
        gr_mat_clear(C, ctxptr);
        gr_mat_clear(D, ctxptr);
        gr_mat_clear(ERR, ctxptr);
        gr_mat_clear(TOL, ctxptr);

        if (ctx == NULL)
            gr_ctx_clear(ctxptr);
    }
}
