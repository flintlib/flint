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

TEST_FUNCTION_START(ca_poly_gcd_euclidean, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B, C, AC, BC, MC, G;
        int success;

        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);
        ca_poly_init(C, ctx);
        ca_poly_init(AC, ctx);
        ca_poly_init(BC, ctx);
        ca_poly_init(MC, ctx);
        ca_poly_init(G, ctx);

        if (n_randint(state, 10) != 0)
        {
            do {
                ca_poly_randtest_rational(A, state, 4, 10, ctx);
                ca_poly_randtest_rational(B, state, 4, 10, ctx);
                success = ca_poly_gcd_euclidean(G, A, B, ctx);
            }
            while (!success || ca_poly_check_is_one(G, ctx) != T_TRUE);

            ca_poly_randtest(C, state, 4, 2, 10, ctx);
        }
        else
        {
            do {
                ca_poly_randtest(A, state, 3, 2, 10, ctx);
                ca_poly_randtest(B, state, 3, 2, 10, ctx);
                success = ca_poly_gcd_euclidean(G, A, B, ctx);
            }
            while (!success || ca_poly_check_is_one(G, ctx) != T_TRUE);

            ca_poly_randtest_rational(C, state, 3, 10, ctx);
        }

        ca_poly_mul(AC, A, C, ctx);
        ca_poly_mul(BC, B, C, ctx);

        switch (n_randint(state, 3))
        {
            case 0:
                ca_poly_set(G, AC, ctx);
                success = ca_poly_gcd_euclidean(G, G, BC, ctx);
                break;
            case 1:
                ca_poly_set(G, BC, ctx);
                success = ca_poly_gcd_euclidean(G, AC, G, ctx);
                break;
            default:
                success = ca_poly_gcd_euclidean(G, AC, BC, ctx);
                break;
        }

        if (success)
        {
            ca_poly_make_monic(MC, C, ctx);

            if (ca_poly_check_equal(MC, G, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                flint_printf("C = "); ca_poly_print(C, ctx); flint_printf("\n");
                flint_printf("G = "); ca_poly_print(G, ctx); flint_printf("\n");

                flint_abort();
            }
        }

        ca_poly_randtest(A, state, 3, 2, 10, ctx);

        success = ca_poly_gcd_euclidean(G, A, A, ctx);

        if (success)
        {
            ca_poly_make_monic(B, A, ctx);

            if (ca_poly_check_equal(G, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL (self)\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                flint_printf("G = "); ca_poly_print(G, ctx); flint_printf("\n");

                flint_abort();
            }
        }

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        ca_poly_clear(C, ctx);
        ca_poly_clear(AC, ctx);
        ca_poly_clear(BC, ctx);
        ca_poly_clear(MC, ctx);
        ca_poly_clear(G, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
