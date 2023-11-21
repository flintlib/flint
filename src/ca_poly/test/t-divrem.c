/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_poly.h"

TEST_FUNCTION_START(ca_poly_divrem, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B, Q, R, T;
        int success;

        /* Test A = Q*B + R */
        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);
        ca_poly_init(Q, ctx);
        ca_poly_init(R, ctx);
        ca_poly_init(T, ctx);

        ca_poly_randtest(A, state, 4, 2, 10, ctx);
        ca_poly_randtest(B, state, 4, 2, 10, ctx);
        ca_poly_randtest(Q, state, 4, 2, 10, ctx);
        ca_poly_randtest(R, state, 4, 2, 10, ctx);
        ca_poly_randtest(T, state, 4, 2, 10, ctx);

        /* test aliasing */
        switch (n_randint(state, 5))
        {
            case 0:
                ca_poly_set(Q, A, ctx);
                success = ca_poly_divrem(Q, R, Q, B, ctx);
                break;
            case 1:
                ca_poly_set(R, A, ctx);
                success = ca_poly_divrem(Q, R, R, B, ctx);
                break;
            case 2:
                ca_poly_set(Q, B, ctx);
                success = ca_poly_divrem(Q, R, A, Q, ctx);
                break;
            case 3:
                ca_poly_set(R, B, ctx);
                success = ca_poly_divrem(Q, R, A, R, ctx);
                break;
            default:
                success = ca_poly_divrem(Q, R, A, B, ctx);
                break;
        }

        if (success)
        {
            ca_poly_mul(T, Q, B, ctx);
            ca_poly_add(T, T, R, ctx);

            if (ca_poly_check_equal(A, T, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                flint_printf("Q = "); ca_poly_print(Q, ctx); flint_printf("\n");
                flint_printf("R = "); ca_poly_print(R, ctx); flint_printf("\n");
                flint_printf("T = "); ca_poly_print(T, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        ca_poly_clear(Q, ctx);
        ca_poly_clear(R, ctx);
        ca_poly_clear(T, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
