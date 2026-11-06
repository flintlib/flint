/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

TEST_GR_FUNCTION_START(gr_mat_conjugate_transpose, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        slong r, c;
        gr_ctx_t ctx;
        gr_mat_t A, B, C, D;
        int alias = n_randint(state, 2);

        gr_ctx_init_random(ctx, state);

        r = n_randint(state, 6);
        c = alias ? r : n_randint(state, 6);

        gr_mat_init(A, r, c, ctx);
        gr_mat_init(B, c, r, ctx);
        gr_mat_init(C, r, c, ctx);
        gr_mat_init(D, c, r, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(B, state, ctx));

        /* B = A^H, C = B^H; expect C == A */
        if (alias)
        {
            status |= gr_mat_set(B, A, ctx);
            status |= gr_mat_conjugate_transpose(B, B, ctx);
        }
        else
        {
            status |= gr_mat_conjugate_transpose(B, A, ctx);
        }

        status |= gr_mat_conjugate_transpose(C, B, ctx);

        if (status == GR_SUCCESS && gr_mat_equal(C, A, ctx) == T_FALSE)
        {
            flint_printf("FAIL: (A^H)^H != A\n");
            gr_ctx_println(ctx);
            flint_printf("alias = %d\n", alias);
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "), gr_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "), gr_mat_print(C, ctx); flint_printf("\n");
            flint_abort();
        }

        /* A^H = conj(A^T) entrywise */
        if (status == GR_SUCCESS)
        {
            slong i, j;
            slong sz = ctx->sizeof_elem;

            status |= gr_mat_transpose(D, A, ctx);
            for (i = 0; i < c; i++)
                for (j = 0; j < r; j++)
                    status |= gr_conj(GR_MAT_ENTRY(D, i, j, sz), GR_MAT_ENTRY(D, i, j, sz), ctx);

            if (status == GR_SUCCESS && gr_mat_equal(D, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL: A^H != conj(A^T)\n");
                gr_ctx_println(ctx);
                flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B = "), gr_mat_print(B, ctx); flint_printf("\n");
                flint_printf("D = "), gr_mat_print(D, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
