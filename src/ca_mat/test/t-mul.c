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

TEST_FUNCTION_START(ca_mat_mul, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, C, AB, BC, X, Y;
        slong M, N, K, L;

        /* Test (A*B)*C = A*(B*C) */
        M = n_randint(state, 4);
        N = n_randint(state, 4);
        K = n_randint(state, 4);
        L = n_randint(state, 4);

        ca_ctx_init(ctx);

        ca_mat_init(A, M, N, ctx);
        ca_mat_init(B, N, K, ctx);
        ca_mat_init(C, K, L, ctx);
        ca_mat_init(AB, M, K, ctx);
        ca_mat_init(BC, N, L, ctx);
        ca_mat_init(X, M, L, ctx);
        ca_mat_init(Y, M, L, ctx);

        ca_mat_randtest(A, state, 2, 10, ctx);
        ca_mat_randtest(B, state, 2, 10, ctx);
        ca_mat_randtest(C, state, 2, 10, ctx);
        ca_mat_randtest(AB, state, 2, 10, ctx);
        ca_mat_randtest(BC, state, 2, 10, ctx);
        ca_mat_randtest(X, state, 2, 10, ctx);
        ca_mat_randtest(Y, state, 2, 10, ctx);

        ca_mat_mul(AB, A, B, ctx);
        ca_mat_mul(X, AB, C, ctx);
        ca_mat_mul(BC, B, C, ctx);
        ca_mat_mul(Y, A, BC, ctx);

        if (ca_mat_check_equal(X, Y, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); ca_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); ca_mat_print(C, ctx); flint_printf("\n");
            flint_printf("X = "); ca_mat_print(X, ctx); flint_printf("\n");
            flint_printf("Y = "); ca_mat_print(Y, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(C, ctx);
        ca_mat_clear(AB, ctx);
        ca_mat_clear(BC, ctx);
        ca_mat_clear(X, ctx);
        ca_mat_clear(Y, ctx);

        ca_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, C, AB, BC, X, Y;
        slong M, N, K, L;

        /* Test (A*B)*C = A*(B*C), rational matrices */
        M = n_randint(state, 5);
        N = n_randint(state, 5);
        K = n_randint(state, 5);
        L = n_randint(state, 5);

        ca_ctx_init(ctx);

        ca_mat_init(A, M, N, ctx);
        ca_mat_init(B, N, K, ctx);
        ca_mat_init(C, K, L, ctx);
        ca_mat_init(AB, M, K, ctx);
        ca_mat_init(BC, N, L, ctx);
        ca_mat_init(X, M, L, ctx);
        ca_mat_init(Y, M, L, ctx);

        ca_mat_randtest_rational(A, state, 200, ctx);
        ca_mat_randtest_rational(B, state, 200, ctx);
        ca_mat_randtest_rational(C, state, 200, ctx);
        ca_mat_randtest_rational(AB, state, 200, ctx);
        ca_mat_randtest_rational(BC, state, 200, ctx);
        ca_mat_randtest_rational(X, state, 200, ctx);
        ca_mat_randtest_rational(Y, state, 200, ctx);

        ca_mat_mul(AB, A, B, ctx);
        ca_mat_mul(X, AB, C, ctx);
        ca_mat_mul(BC, B, C, ctx);
        ca_mat_mul(Y, A, BC, ctx);

        if (ca_mat_check_equal(X, Y, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); ca_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); ca_mat_print(C, ctx); flint_printf("\n");
            flint_printf("X = "); ca_mat_print(X, ctx); flint_printf("\n");
            flint_printf("Y = "); ca_mat_print(Y, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(C, ctx);
        ca_mat_clear(AB, ctx);
        ca_mat_clear(BC, ctx);
        ca_mat_clear(X, ctx);
        ca_mat_clear(Y, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
