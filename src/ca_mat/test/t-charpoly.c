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

TEST_FUNCTION_START(ca_mat_charpoly, state)
{
    slong iter;

    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, fA;
        ca_poly_t f;
        slong n;

        ca_ctx_init(ctx);

        n = n_randint(state, 5);
        ca_mat_init(A, n, n, ctx);
        ca_mat_init(fA, n, n, ctx);
        ca_poly_init(f, ctx);

        ca_mat_randtest(A, state, 1, 5, ctx);

        ca_mat_charpoly(f, A, ctx);

        ca_mat_ca_poly_evaluate(fA, f, A, ctx);

        if (f->length != (n + 1) || ca_mat_check_is_zero(fA, ctx) != T_TRUE)
        {
            flint_printf("FAIL\n");
            flint_printf("A = "), ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("f = "), ca_poly_print(f, ctx); flint_printf("\n");
            flint_printf("fA = "); ca_mat_print(fA, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(fA, ctx);
        ca_poly_clear(f, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
