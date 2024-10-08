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

void gr_mat_test_det(gr_method_mat_unary_op_get_scalar det_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t ctx2;
        gr_ctx_ptr ctxptr;

        if (ctx == NULL)
        {
            gr_ctx_init_random(ctx2, state);
            ctxptr = ctx2;
        }
        else
            ctxptr = ctx;

        {
            gr_mat_t A, B, AB;
            gr_ptr detA, detB, detAB, detAdetB;
            slong n;
            int status = GR_SUCCESS;

            n = n_randint(state, maxn + 1);

            gr_mat_init(A, n, n, ctxptr);
            gr_mat_init(B, n, n, ctxptr);
            gr_mat_init(AB, n, n, ctxptr);

            detA = gr_heap_init(ctxptr);
            detB = gr_heap_init(ctxptr);
            detAB = gr_heap_init(ctxptr);
            detAdetB = gr_heap_init(ctxptr);

            status |= gr_mat_randtest(A, state, ctxptr);
            status |= gr_mat_randtest(B, state, ctxptr);
            status |= gr_mat_mul(AB, A, B, ctxptr);

            status |= det_impl(detA, A, ctxptr);
            status |= det_impl(detB, B, ctxptr);
            status |= det_impl(detAB, AB, ctxptr);
            status |= gr_mul(detAdetB, detA, detB, ctxptr);

            /* Check that the output isn't just 0 */
            if (status == GR_SUCCESS && gr_mat_is_one(A, ctxptr) == T_TRUE &&
                  gr_is_one(detA, ctxptr) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                gr_ctx_println(ctxptr);
                flint_printf("A = "); gr_mat_print(A, ctxptr); flint_printf("\n");
                flint_printf("detA = "); gr_print(detA, ctxptr); flint_printf("\n");
                flint_abort();
            }

            if (status == GR_SUCCESS && gr_equal(detAB, detAdetB, ctxptr) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                gr_ctx_println(ctxptr);
                flint_printf("A = "); gr_mat_print(A, ctxptr); flint_printf("\n");
                flint_printf("B = "); gr_mat_print(B, ctxptr); flint_printf("\n");
                flint_printf("AB = "); gr_mat_print(AB, ctxptr); flint_printf("\n");
                flint_printf("detA = "); gr_print(detA, ctxptr); flint_printf("\n");
                flint_printf("detB = "); gr_print(detB, ctxptr); flint_printf("\n");
                flint_printf("detAB = "); gr_print(detAB, ctxptr); flint_printf("\n");
                flint_printf("detAdetB = "); gr_print(detAdetB, ctxptr); flint_printf("\n");
                flint_abort();
            }

            if ((status & GR_DOMAIN) && !(status & GR_UNABLE))
            {
                flint_printf("FAIL (flags)\n\n");
                gr_ctx_println(ctxptr);
                flint_printf("A = "); gr_mat_print(A, ctxptr); flint_printf("\n");
                flint_printf("B = "); gr_mat_print(B, ctxptr); flint_printf("\n");
                flint_printf("AB = "); gr_mat_print(AB, ctxptr); flint_printf("\n");
                flint_abort();
            }

            gr_mat_clear(A, ctxptr);
            gr_mat_clear(B, ctxptr);
            gr_mat_clear(AB, ctxptr);

            /* Check error handling */
            gr_mat_init(A, n, n + 1, ctxptr);
            status = det_impl(detA, A, ctxptr);

            if (status == GR_SUCCESS)
            {
                flint_printf("FAIL (nonsquare matrix)\n\n");
                gr_ctx_println(ctxptr);
                flint_printf("A = "); gr_mat_print(A, ctxptr); flint_printf("\n");
                flint_abort();
            }

            gr_mat_clear(A, ctxptr);

            gr_heap_clear(detA, ctxptr);
            gr_heap_clear(detB, ctxptr);
            gr_heap_clear(detAB, ctxptr);
            gr_heap_clear(detAdetB, ctxptr);
        }

        if (ctx == NULL)
            gr_ctx_clear(ctxptr);
    }
}
