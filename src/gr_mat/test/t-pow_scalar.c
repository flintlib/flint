/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_GR_FUNCTION_START(gr_mat_pow_scalar, state, count_success, count_domain, count_unable)
{
    slong iter;

    /* Check A^c * A = A^(c+1) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong n;
        gr_ctx_t ctx;
        gr_ptr c, c1;
        gr_mat_t A, Ac, Ac1, AcA;
        int status = GR_SUCCESS;
        int which;

        while (1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_field(ctx) == T_TRUE)
                break;
            gr_ctx_clear(ctx);
        }

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 3);
        else
            n = n_randint(state, 6);

        which = n_randint(state, 2);

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(Ac, n, n, ctx);
        gr_mat_init(Ac1, n, n, ctx);
        gr_mat_init(AcA, n, n, ctx);
        c = gr_heap_init(ctx);
        c1 = gr_heap_init(ctx);

        if (n_randint(state, 4) == 0)
        {
            GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        }
        else
        {
            fmpq_mat_t Q;
            fmpq_mat_init(Q, n, n);
            fmpq_mat_randtest(Q, state, 5);
            status = gr_mat_set_fmpq_mat(A, Q, ctx);
            fmpq_mat_clear(Q);
        }

        GR_MUST_SUCCEED(gr_randtest_small(c, state, ctx));

        if (which)
            status = gr_mat_pow_scalar(Ac, A, c, ctx);
        else
            status = gr_mat_pow_scalar_jordan(Ac, A, c, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_add_ui(c1, c, 1, ctx);
            status |= gr_mat_pow_scalar(Ac1, A, c1, ctx);
            status |= gr_mat_mul(AcA, Ac, A, ctx);

            if (status == GR_SUCCESS && gr_mat_equal(AcA, Ac1, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n");
                gr_ctx_println(ctx);
                flint_printf("c = "), gr_print(c, ctx); flint_printf("\n");
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("Ac = "), gr_mat_print(Ac, ctx); flint_printf("\n");
                flint_printf("Ac1 = "), gr_mat_print(Ac1, ctx); flint_printf("\n");
                flint_printf("AcA = "), gr_mat_print(AcA, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(Ac, ctx);
        gr_mat_clear(Ac1, ctx);
        gr_mat_clear(AcA, ctx);
        gr_heap_clear(c, ctx);
        gr_heap_clear(c1, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
