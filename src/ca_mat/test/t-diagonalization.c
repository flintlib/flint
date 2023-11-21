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

TEST_FUNCTION_START(ca_mat_diagonalization, state)
{
    slong iter;

    /* todo: test non-diagonalizable matrices */
    /* todo: test aliasing */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, P, D, Pinv, B;
        slong n;
        truth_t result;

        ca_ctx_init(ctx);
        n = n_randint(state, 4);

        ca_mat_init(A, n, n, ctx);
        ca_mat_init(P, n, n, ctx);
        ca_mat_init(D, n, n, ctx);
        ca_mat_init(Pinv, n, n, ctx);
        ca_mat_init(B, n, n, ctx);

        if (n <= 2)
            ca_mat_randtest(A, state, 1, 5, ctx);
        else
            ca_mat_randtest_rational(A, state, 5, ctx);

        ca_mat_randtest(D, state, 1, 5, ctx);
        ca_mat_randtest(P, state, 1, 5, ctx);

        result = ca_mat_diagonalization(D, P, A, ctx);

        if (result == T_TRUE)
        {
            result = ca_mat_inv(Pinv, P, ctx);

            if (result == T_FALSE)
            {
                flint_printf("FAIL (inversion of P)\n");
                flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("D: "); ca_mat_print(D, ctx); flint_printf("\n");
                flint_printf("P: "); ca_mat_print(P, ctx); flint_printf("\n");
                flint_abort();
            }

            if (result == T_TRUE)
            {
                ca_mat_mul(B, P, D, ctx);
                ca_mat_mul(B, B, Pinv, ctx);

                result = ca_mat_check_equal(A, B, ctx);

                if (result == T_FALSE)
                {
                    flint_printf("FAIL\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("D: "); ca_mat_print(D, ctx); flint_printf("\n");
                    flint_printf("P: "); ca_mat_print(P, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(D, ctx);
        ca_mat_clear(P, ctx);
        ca_mat_clear(Pinv, ctx);
        ca_mat_clear(B, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
