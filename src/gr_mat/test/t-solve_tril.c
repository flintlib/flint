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

TEST_GR_FUNCTION_START(gr_mat_solve_tril, state, count_success, count_unable, count_domain)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, X, B, Y;
        slong rows, cols, i, j;
        int unit;
        int status = GR_SUCCESS;
        slong sz;

        rows = n_randint(state, 15);
        cols = n_randint(state, 15);
        unit = n_randint(state, 2);

        /* Hack: avoid because slow */
        gr_ctx_init_random(ctx, state);
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        sz = ctx->sizeof_elem;

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
            status |= gr_mat_nonsingular_solve_tril(Y, A, B, unit, ctx);
        }
        else
        {
            status |= gr_mat_set(Y, B, ctx);
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

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(X, ctx);
        gr_mat_clear(Y, ctx);
        gr_ctx_clear(ctx);

    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
