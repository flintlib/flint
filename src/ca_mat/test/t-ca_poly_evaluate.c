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

TEST_FUNCTION_START(ca_mat_ca_poly_evaluate, state)
{
    slong iter;

    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, fA, gA, fgA, fAgA;
        ca_poly_t f, g, fg;
        slong n;

        n = n_randint(state, 3);

        ca_ctx_init(ctx);

        ca_mat_init(A, n, n, ctx);
        ca_mat_init(fA, n, n, ctx);
        ca_mat_init(gA, n, n, ctx);
        ca_mat_init(fgA, n, n, ctx);
        ca_mat_init(fAgA, n, n, ctx);

        ca_poly_init(f, ctx);
        ca_poly_init(g, ctx);
        ca_poly_init(fg, ctx);

        ca_mat_randtest(A, state, 1, 10, ctx);
        ca_mat_randtest(fA, state, 2, 10, ctx);
        ca_mat_randtest(gA, state, 2, 10, ctx);
        ca_mat_randtest(fgA, state, 2, 10, ctx);

        ca_poly_randtest_rational(f, state, 10, 5, ctx);
        ca_poly_randtest_rational(g, state, 10, 5, ctx);
        ca_poly_add(fg, f, g, ctx);

        ca_mat_ca_poly_evaluate(fA, f, A, ctx);
        ca_mat_set(gA, A, ctx);
        ca_mat_ca_poly_evaluate(gA, g, gA, ctx);

        ca_mat_ca_poly_evaluate(fgA, fg, A, ctx);

        ca_mat_add(fAgA, fA, gA, ctx);

        if (ca_mat_check_equal(fgA, fAgA, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("f = "); ca_poly_print(f, ctx); flint_printf("\n");
            flint_printf("g = "); ca_poly_print(g, ctx); flint_printf("\n");
            flint_printf("fgA = "); ca_mat_print(fgA, ctx); flint_printf("\n");
            flint_printf("fAgA = "); ca_mat_print(fAgA, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(fA, ctx);
        ca_mat_clear(gA, ctx);
        ca_mat_clear(fgA, ctx);
        ca_mat_clear(fAgA, ctx);

        ca_poly_clear(f, ctx);
        ca_poly_clear(g, ctx);
        ca_poly_clear(fg, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
