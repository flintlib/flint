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

TEST_FUNCTION_START(ca_poly_roots, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_vec_t R;
        ca_poly_t A, B;
        ulong * exp;

        ca_ctx_init(ctx);

        ca_vec_init(R, 0, ctx);
        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);

        if (n_randint(state, 10) != 0)
        {
            ca_poly_randtest(A, state, 5, 1, 10, ctx);
        }
        else
        {
            ca_poly_randtest(A, state, 3, 0, 5, ctx);
            ca_poly_randtest_rational(B, state, 3, 5, ctx);
            ca_poly_mul(A, A, B, ctx);
            if (n_randint(state, 2))
                ca_poly_mul(A, A, B, ctx);
            if (n_randint(state, 2))
                ca_poly_mul(A, A, B, ctx);
        }

        exp = flint_malloc(sizeof(ulong) * A->length);

        if (ca_poly_roots(R, exp, A, ctx))
        {
            ca_poly_set_roots(B, R, exp, ctx);

            if (A->length)
                ca_poly_mul_ca(B, B, A->coeffs + A->length - 1, ctx);

            if (ca_poly_check_equal(A, B, ctx) == T_FALSE)
            {
                slong i;
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                flint_printf("R = "); ca_vec_print(R, ctx); flint_printf("\n");
                for (i = 0; i < R->length; i++)
                    flint_printf("e[%wd] = %wu\n", i, exp[i]);
                flint_printf("\n");
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                flint_abort();
            }

            if (0)
            {
                printf("=================================================================\n\n");
                printf("EQUAL: "); truth_print(ca_poly_check_equal(A, B, ctx)); printf("\n\n");
                flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n\n");
                flint_printf("R = "); ca_vec_print(R, ctx); flint_printf("\n\n");
                ca_poly_sub(B, A, B, ctx);
                flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n\n");
            }
        }

        ca_vec_clear(R, ctx);
        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        flint_free(exp);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
