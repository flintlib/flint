/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_mat.h"

TEST_FUNCTION_START(ca_mat_solve_triu, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, X, B, Y;
        slong rows, cols, i, j;
        int unit;

        rows = n_randint(state, 15);
        cols = n_randint(state, 15);
        unit = n_randint(state, 2);

        ca_ctx_init(ctx);
        ca_mat_init(A, rows, rows, ctx);
        ca_mat_init(B, rows, cols, ctx);
        ca_mat_init(X, rows, cols, ctx);
        ca_mat_init(Y, rows, cols, ctx);

        /* todo: test number fields, etc. */
        ca_mat_randtest_rational(A, state, 5, ctx);
        ca_mat_randtest_rational(X, state, 5, ctx);
        ca_mat_randtest_rational(Y, state, 5, ctx);

        for (i = 0; i < rows; i++)
        {
            if (unit)
                ca_one(ca_mat_entry(A, i, i), ctx);
            else
                ca_set_ui(ca_mat_entry(A, i, i), 1 + n_randint(state, 100), ctx);

            for (j = 0; j < i; j++)
                ca_zero(ca_mat_entry(A, i, j), ctx);
        }

        ca_mat_mul(B, A, X, ctx);

        if (unit)  /* check that diagonal entries are ignored */
        {
            for (i = 0; i < rows; i++)
                ca_set_ui(ca_mat_entry(A, i, i), 1 + n_randint(state, 100), ctx);
        }

        /* Check Y = A^(-1) * (A * X) = X */
        if (n_randint(state, 2))
        {
            ca_mat_solve_triu(Y, A, B, unit, ctx);
        }
        else
        {
            ca_mat_set(Y, B, ctx);
            ca_mat_solve_triu(Y, A, Y, unit, ctx);
        }

        if (ca_mat_check_equal(Y, X, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            flint_printf("A = \n"); ca_mat_print(A, ctx); flint_printf("\n\n");
            flint_printf("B = \n"); ca_mat_print(B, ctx); flint_printf("\n\n");
            flint_printf("X = \n"); ca_mat_print(X, ctx); flint_printf("\n\n");
            flint_printf("Y = \n"); ca_mat_print(Y, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(X, ctx);
        ca_mat_clear(Y, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
