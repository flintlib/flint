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

TEST_FUNCTION_START(ca_mat_companion, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A;
        ca_poly_t f, g;
        slong n;

        ca_ctx_init(ctx);

        ca_poly_init(f, ctx);
        ca_poly_init(g, ctx);

        do {
            ca_poly_randtest(f, state, 1 + n_randint(state, 4), 2, 5, ctx);
        } while (f->length < 1);

        n = f->length - 1;

        ca_mat_init(A, n, n, ctx);
        ca_mat_randtest_rational(A, state, 100, ctx);

        if (ca_mat_companion(A, f, ctx))
        {
            ca_mat_charpoly(g, A, ctx);
            ca_poly_mul_ca(g, g, f->coeffs + n, ctx);

            if (ca_poly_check_equal(g, f, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n");
                flint_printf("A = "), ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("f = "), ca_poly_print(f, ctx); flint_printf("\n");
                flint_printf("g = "); ca_poly_print(g, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        ca_mat_clear(A, ctx);
        ca_poly_clear(f, ctx);
        ca_poly_clear(g, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
