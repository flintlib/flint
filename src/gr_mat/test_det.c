/*
    Copyright (C) 2022, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you grn redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

void gr_mat_test_det(gr_method_mat_unary_op_get_scalar det_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
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
            gr_mat_t A, B, AB;
            gr_ptr detA, detB, detAB, detAdetB;
            slong n;
            int status = GR_SUCCESS;

            n = n_randint(state, maxn + 1);

            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(AB, n, n, ctx);

            detA = gr_heap_init(ctx);
            detB = gr_heap_init(ctx);
            detAB = gr_heap_init(ctx);
            detAdetB = gr_heap_init(ctx);

            status |= gr_mat_randtest(A, state, ctx);
            status |= gr_mat_randtest(B, state, ctx);
            status |= gr_mat_mul(AB, A, B, ctx);

            status |= det_impl(detA, A, ctx);
            status |= det_impl(detB, B, ctx);
            status |= det_impl(detAB, AB, ctx);
            status |= gr_mul(detAdetB, detA, detB, ctx);

            /* Check that the output isn't just 0 */
            if (status == GR_SUCCESS && gr_mat_is_one(A, ctx) == T_TRUE &&
                  gr_is_one(detA, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("detA = "); gr_print(detA, ctx); flint_printf("\n");
                flint_abort();
            }

            if (status == GR_SUCCESS && gr_equal(detAB, detAdetB, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); gr_mat_print(B, ctx); flint_printf("\n");
                flint_printf("AB = "); gr_mat_print(AB, ctx); flint_printf("\n");
                flint_printf("detA = "); gr_print(detA, ctx); flint_printf("\n");
                flint_printf("detB = "); gr_print(detB, ctx); flint_printf("\n");
                flint_printf("detAB = "); gr_print(detAB, ctx); flint_printf("\n");
                flint_printf("detAdetB = "); gr_print(detAdetB, ctx); flint_printf("\n");
                flint_abort();
            }

            if ((status & GR_DOMAIN) && !(status & GR_UNABLE))
            {
                flint_printf("FAIL (flags)\n\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); gr_mat_print(B, ctx); flint_printf("\n");
                flint_printf("AB = "); gr_mat_print(AB, ctx); flint_printf("\n");
                flint_abort();
            }

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(AB, ctx);

            /* Check error handling */
            gr_mat_init(A, n, n + 1, ctx);
            status = det_impl(detA, A, ctx);

            if (status == GR_SUCCESS)
            {
                flint_printf("FAIL (nonsquare matrix)\n\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_abort();
            }

            gr_mat_clear(A, ctx);

            gr_heap_clear(detA, ctx);
            gr_heap_clear(detB, ctx);
            gr_heap_clear(detAB, ctx);
            gr_heap_clear(detAdetB, ctx);
        }

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}
