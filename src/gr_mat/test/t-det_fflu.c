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
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_GR_FUNCTION_START(gr_mat_det_fflu, state, count_success, count_unable, count_domain)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status = GR_SUCCESS;
        slong n;
        gr_ctx_t ctx;
        gr_mat_t A, B, AB;
        gr_ptr detA, detB, detAB, detAdetB;

        /* Hack: avoid because slow */
        gr_ctx_init_random(ctx, state);
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        n = n_randint(state, 7);

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

        status |= gr_mat_det_fflu(detA, A, ctx);
        status |= gr_mat_det_fflu(detB, B, ctx);
        status |= gr_mat_det_fflu(detAB, AB, ctx);
        status |= gr_mul(detAdetB, detA, detB, ctx);

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

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(AB, ctx);

        gr_heap_clear(detA, ctx);
        gr_heap_clear(detB, ctx);
        gr_heap_clear(detAB, ctx);
        gr_heap_clear(detAdetB, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
