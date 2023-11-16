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

TEST_GR_FUNCTION_START(gr_mat_adjugate, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status = GR_SUCCESS;
        slong n;
        gr_ctx_t ctx;
        gr_mat_t A, B, C;
        gr_ptr d, e;

        gr_ctx_init_random(ctx, state);

        /* Hack: avoid because slow */
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        n = n_randint(state, 5);
        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C, n, n, ctx);

        d = gr_heap_init(ctx);
        e = gr_heap_init(ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

        status |= gr_mat_adjugate_cofactor(B, d, A, ctx);
        status |= gr_mat_adjugate_charpoly(C, e, A, ctx);

        if (status == GR_SUCCESS && (gr_mat_equal(B, C, ctx) == T_FALSE || gr_equal(d, e, ctx) == T_FALSE))
        {
            flint_printf("FAIL\n");
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "), gr_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "), gr_mat_print(C, ctx); flint_printf("\n");
            flint_printf("d = "), gr_print(d, ctx); flint_printf("\n");
            flint_printf("e = "), gr_print(e, ctx); flint_printf("\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);

        gr_heap_clear(d, ctx);
        gr_heap_clear(e, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
