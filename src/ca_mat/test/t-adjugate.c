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

TEST_FUNCTION_START(ca_mat_adjugate, state)
{
    slong iter;

    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, C;
        ca_t d, e;
        slong n;

        ca_ctx_init(ctx);

        n = n_randint(state, 5);
        ca_mat_init(A, n, n, ctx);
        ca_mat_init(B, n, n, ctx);
        ca_mat_init(C, n, n, ctx);
        ca_init(d, ctx);
        ca_init(e, ctx);

        ca_mat_randtest(A, state, 1, 5, ctx);

        ca_mat_adjugate_cofactor(B, d, A, ctx);
        ca_mat_adjugate_charpoly(C, e, A, ctx);

        if (ca_mat_check_equal(B, C, ctx) == T_FALSE || ca_check_equal(d, e, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            flint_printf("A = "), ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "), ca_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "), ca_mat_print(C, ctx); flint_printf("\n");
            flint_printf("d = "), ca_print(d, ctx); flint_printf("\n");
            flint_printf("e = "), ca_print(e, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(C, ctx);
        ca_clear(d, ctx);
        ca_clear(e, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
