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

TEST_FUNCTION_START(ca_mat_dft, state)
{
    slong iter;

    for (iter = 0; iter < 1; iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, C;
        slong n;

        ca_ctx_init(ctx);

        for (n = 0; n <= 16; n++)
        {
            ca_mat_init(A, n, n, ctx);
            ca_mat_init(B, n, n, ctx);
            ca_mat_init(C, n, n, ctx);

            ca_mat_dft(A, 0, ctx);
            ca_mat_dft(B, 1, ctx);
            ca_mat_mul(C, A, B, ctx);

            if (ca_mat_check_is_one(C, ctx) != T_TRUE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_mat_print(B, ctx); flint_printf("\n");
                flint_printf("C = "); ca_mat_print(C, ctx); flint_printf("\n");
                flint_abort();
            }

            ca_mat_dft(A, 2, ctx);
            ca_mat_dft(B, 3, ctx);
            ca_mat_mul(C, A, B, ctx);

            if (ca_mat_check_is_one(C, ctx) != T_TRUE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); ca_mat_print(B, ctx); flint_printf("\n");
                flint_printf("C = "); ca_mat_print(C, ctx); flint_printf("\n");
                flint_abort();
            }

            ca_mat_clear(A, ctx);
            ca_mat_clear(B, ctx);
            ca_mat_clear(C, ctx);
        }

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
