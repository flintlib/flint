/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "gr_vec.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_diagonalization, state)
{
    slong iter;

    /* todo: test non-diagonalizable matrices */
    /* todo: test aliasing */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, L, R, B, LR;
        gr_vec_t D;
        slong n;
        int status = GR_SUCCESS;
        int haveL, haveR;

        gr_ctx_init_random(ctx, state);
        n = n_randint(state, 5);
        haveL = n_randint(state, 2);
        haveR = n_randint(state, 2);

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(L, n, n, ctx);
        gr_mat_init(R, n, n, ctx);
        gr_vec_init(D, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(LR, n, n, ctx);

        if (n <= 2)
        {
            status |= gr_mat_randtest(A, state, ctx);
        }
        else
        {
            fmpq_mat_t T;
            fmpq_mat_init(T, n, n);
            fmpq_mat_randtest(T, state, 1);
            status |= gr_mat_set_fmpq_mat(A, T, ctx);
            fmpq_mat_clear(T);
        }
        /* status |= gr_vec_randtest(D, state, ctx); */
        status |= gr_mat_randtest(L, state, ctx);
        status |= gr_mat_randtest(R, state, ctx);

        status |= gr_mat_diagonalization(D, haveL ? L : NULL, haveR ? R : NULL, A, 0, ctx);

        if (status == GR_SUCCESS)
        {
            if (haveL && haveR)
            {
                status = gr_mat_mul(LR, L, R, ctx);

                if (status == GR_SUCCESS && gr_mat_is_one(LR, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: L*R != 1\n");
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("D: "); gr_vec_print(D, ctx); flint_printf("\n");
                    flint_printf("L: "); gr_mat_print(L, ctx); flint_printf("\n");
                    flint_printf("R: "); gr_mat_print(R, ctx); flint_printf("\n");
                    flint_printf("LR: "); gr_mat_print(LR, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
            else if (haveL || haveR)
            {
                if (haveL)
                    status |= gr_mat_inv(R, L, ctx);
                else
                    status |= gr_mat_inv(L, R, ctx);

                if (status == GR_DOMAIN && gr_ctx_is_field(ctx) == T_TRUE)
                {
                    flint_printf("FAIL (inversion of R or L)\n");
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("D: "); gr_vec_print(D, ctx); flint_printf("\n");
                    flint_printf("L: "); gr_mat_print(L, ctx); flint_printf("\n");
                    flint_printf("R: "); gr_mat_print(R, ctx); flint_printf("\n");
                    flint_abort();
                }
            }

            if ((haveL || haveR) && status == GR_SUCCESS)
            {
                /*
                gr_ctx_println(ctx);
                gr_mat_print(A, ctx); printf("\n\n");
                */

                status |= gr_mat_mul_diag(B, R, D, ctx);
                status |= gr_mat_mul(B, B, L, ctx);

                if (status == GR_SUCCESS && gr_mat_equal(A, B, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: RDL != A\n");
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("D: "); gr_vec_print(D, ctx); flint_printf("\n");
                    flint_printf("L: "); gr_mat_print(L, ctx); flint_printf("\n");
                    flint_printf("R: "); gr_mat_print(R, ctx); flint_printf("\n");
                    flint_printf("B: "); gr_mat_print(B, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        gr_mat_clear(A, ctx);
        gr_vec_clear(D, ctx);
        gr_mat_clear(L, ctx);
        gr_mat_clear(R, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(LR, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
