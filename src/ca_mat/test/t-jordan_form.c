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

TEST_FUNCTION_START(ca_mat_jordan_form, state)
{
    slong iter;

    /* Test P J P^{-1} = A for random matrices A generated from Jordan blocks */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, P, Q, B, J;
        slong n, i, j, num_lambda;
        int success;
        slong num_blocks1, offset;
        ca_vec_t lambda1;
        slong *block_lambda1, *block_size1;
        truth_t inv;

        ca_ctx_init(ctx);

        num_lambda = 1 + n_randint(state, 3);
        ca_vec_init(lambda1, num_lambda, ctx);

        offset = n_randint(state, 2);
        for (i = 0; i < num_lambda; i++)
            ca_set_ui(ca_vec_entry(lambda1, i), offset + i, ctx);

        block_lambda1 = flint_malloc(sizeof(slong) * 4 * num_lambda);
        block_size1 = flint_malloc(sizeof(slong) * 4 * num_lambda);

        num_blocks1 = 0;
        n = 0;
        for (i = 0; i < num_lambda; i++)
        {
            slong blocks = 1 + n_randint(state, 3);

            for (j = 0; j < blocks; j++)
            {
                block_size1[num_blocks1] = 1 + n_randint(state, 3);
                block_lambda1[num_blocks1] = i;
                n += block_size1[num_blocks1];
                num_blocks1++;
            }
        }

        ca_mat_init(A, n, n, ctx);
        ca_mat_init(P, n, n, ctx);
        ca_mat_init(Q, n, n, ctx);
        ca_mat_init(B, n, n, ctx);
        ca_mat_init(J, n, n, ctx);

        ca_mat_set_jordan_blocks(A, lambda1, num_blocks1, block_lambda1, block_size1, ctx);

        do {
            ca_mat_randtest_rational(P, state, 1, ctx);
        } while (ca_mat_inv(Q, P, ctx) != T_TRUE);

        ca_mat_mul(B, P, A, ctx);
        ca_mat_mul(B, B, Q, ctx);

        ca_mat_zero(P, ctx);
        ca_mat_zero(Q, ctx);

        success = ca_mat_jordan_form(J, P, A, ctx);

        if (success)
        {
            inv = ca_mat_inv(Q, P, ctx);

            if (inv == T_FALSE)
            {
                flint_printf("FAIL (matrix not invertible)\n");
                flint_printf("A: "); ca_mat_printn(A, 5, ctx); flint_printf("\n");
                flint_printf("J: "); ca_mat_printn(J, 5, ctx); flint_printf("\n");
                flint_printf("P: "); ca_mat_printn(P, 5, ctx); flint_printf("\n");
                flint_abort();
            }

            ca_mat_mul(B, P, J, ctx);
            ca_mat_mul(B, B, Q, ctx);

            if (ca_mat_check_equal(B, A, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n");
                flint_printf("A: "); ca_mat_printn(A, 5, ctx); flint_printf("\n");
                flint_printf("B: "); ca_mat_printn(B, 5, ctx); flint_printf("\n");
                flint_printf("J: "); ca_mat_printn(J, 5, ctx); flint_printf("\n");
                flint_printf("P: "); ca_mat_printn(P, 5, ctx); flint_printf("\n");
                flint_printf("Q: "); ca_mat_printn(Q, 5, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(P, ctx);
        ca_mat_clear(Q, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(J, ctx);
        ca_vec_clear(lambda1, ctx);
        flint_free(block_lambda1);
        flint_free(block_size1);

        ca_ctx_clear(ctx);
    }

    /* Test P J P^{-1} = A for random matrices A */
    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, J, P, Q;
        slong n;
        int success;
        truth_t inv;

        n = n_randint(state, 5);

        ca_ctx_init(ctx);
        ca_mat_init(A, n, n, ctx);
        ca_mat_init(B, n, n, ctx);
        ca_mat_init(J, n, n, ctx);
        ca_mat_init(P, n, n, ctx);
        ca_mat_init(Q, n, n, ctx);

        if (n <= 2 && n_randint(state, 2) == 0)
            ca_mat_randtest(A, state, 1, 5, ctx);
        else
            ca_mat_randtest_rational(A, state, 5, ctx);

        /* test aliasing */
        if (n_randint(state, 2))
            success = ca_mat_jordan_form(J, P, A, ctx);
        else
        {
            ca_mat_set(P, A, ctx);
            success = ca_mat_jordan_form(J, P, P, ctx);
        }

        if (success)
        {
            inv = ca_mat_inv(Q, P, ctx);

            if (inv == T_FALSE)
            {
                flint_printf("FAIL (matrix not invertible)\n");
                flint_printf("A: "); ca_mat_printn(A, 5, ctx); flint_printf("\n");
                flint_printf("J: "); ca_mat_printn(J, 5, ctx); flint_printf("\n");
                flint_printf("P: "); ca_mat_printn(P, 5, ctx); flint_printf("\n");
                flint_abort();
            }

            ca_mat_mul(B, P, J, ctx);
            ca_mat_mul(B, B, Q, ctx);

            if (ca_mat_check_equal(B, A, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n");
                flint_printf("A: "); ca_mat_printn(A, 5, ctx); flint_printf("\n");
                flint_printf("B: "); ca_mat_printn(B, 5, ctx); flint_printf("\n");
                flint_printf("J: "); ca_mat_printn(J, 5, ctx); flint_printf("\n");
                flint_printf("P: "); ca_mat_printn(P, 5, ctx); flint_printf("\n");
                flint_printf("Q: "); ca_mat_printn(Q, 5, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(J, ctx);
        ca_mat_clear(P, ctx);
        ca_mat_clear(Q, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
