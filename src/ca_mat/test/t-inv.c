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

TEST_FUNCTION_START(ca_mat_inv, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, C;
        slong n;
        truth_t success;

        ca_ctx_init(ctx);
        n = n_randint(state, 6);

        ca_mat_init(A, n, n, ctx);
        ca_mat_init(B, n, n, ctx);
        ca_mat_init(C, n, n, ctx);

        if (n_randint(state, 2) || n > 3)
            ca_mat_randtest_rational(A, state, 5, ctx);
        else
            ca_mat_randtest(A, state, 1, 5, ctx);

        success = ca_mat_inv(B, A, ctx);

        if (success == T_TRUE)
        {
            ca_mat_mul(C, A, B, ctx);

            if (ca_mat_check_is_one(C, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_mat_print(B, ctx); flint_printf("\n");
                flint_printf("C = "); ca_mat_print(C, ctx); flint_printf("\n");
                flint_abort();
            }
        }
        else if (success == T_FALSE)
        {
            ca_t det;
            ca_init(det, ctx);

            ca_mat_det(det, A, ctx);

            if (ca_check_is_zero(det, ctx) == T_FALSE)
            {
                flint_printf("FAIL (singular)\n\n");
                flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("det = "); ca_print(det, ctx); flint_printf("\n");
                flint_abort();
            }

            ca_clear(det, ctx);
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(C, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
