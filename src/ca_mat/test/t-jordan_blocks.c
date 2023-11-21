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

int
check_jordan_forms(ca_vec_t lambda1, slong num_blocks1, slong * block_lambda1, slong * block_size1,
                   ca_vec_t lambda2, slong num_blocks2, slong * block_lambda2, slong * block_size2,
                    ca_ctx_t ctx)
{
    slong i, j;
    slong * used;
    slong found;
    int result;

    if (num_blocks1 != num_blocks2)
        return 0;

    if (ca_vec_length(lambda1, ctx) != ca_vec_length(lambda2, ctx))
        return 0;

    /* note: not ordered
    for (i = 0; i < ca_vec_length(lambda1, ctx); i++)
        if (ca_check_equal(ca_vec_entry(lambda1, i), ca_vec_entry(lambda2, i), ctx) == T_FALSE)
            return 0;
    */

    used = flint_calloc(sizeof(slong), num_blocks1);
    result = 1;

    for (i = 0; result && i < num_blocks1; i++)
    {
        found = -1;

        for (j = 0; found == -1 && j < num_blocks2; j++)
        {
            if (!used[j])
            {
                if (/* block_lambda1[i] == block_lambda2[j] -- if ordered */
                    ca_check_equal(ca_vec_entry(lambda1, block_lambda1[i]), ca_vec_entry(lambda2, block_lambda2[j]), ctx) != T_FALSE
                        && block_size1[i] == block_size2[j])
                {
                    found = j;
                    used[j] = 1;
                }
            }
        }

        if (found == -1)
            result = 0;
    }

    flint_free(used);
    return result;
}

TEST_FUNCTION_START(ca_mat_jordan_blocks, state)
{
    slong iter;

    /* Test J(A) = J(P * A * P^-1) for random rational Jordan block matrix A and random P */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, P, Q, B;
        slong n, i, j, num_lambda;
        int success1, success2;
        slong num_blocks1, num_blocks2, offset;
        ca_vec_t lambda1, lambda2;
        slong *block_lambda1, *block_size1, *block_lambda2, *block_size2;

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
        ca_vec_init(lambda2, 0, ctx);

        block_lambda2 = flint_malloc(sizeof(slong) * n);
        block_size2 = flint_malloc(sizeof(slong) * n);

        ca_mat_set_jordan_blocks(A, lambda1, num_blocks1, block_lambda1, block_size1, ctx);

        do {
            ca_mat_randtest_rational(P, state, 1, ctx);
        } while (ca_mat_inv(Q, P, ctx) != T_TRUE);

        ca_mat_mul(B, P, A, ctx);
        ca_mat_mul(B, B, Q, ctx);

        success1 = 1;
        success2 = ca_mat_jordan_blocks(lambda2, &num_blocks2, block_lambda2, block_size2, B, ctx);

        if (success1 && success2)
        {
            if (!check_jordan_forms(lambda1, num_blocks1, block_lambda1, block_size1,
                                    lambda2, num_blocks2, block_lambda2, block_size2, ctx))
            {
                flint_printf("FAIL (different Jordan forms)\n");
                flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B: "); ca_mat_print(B, ctx); flint_printf("\n");
                flint_printf("lambda1: "); ca_vec_print(lambda1, ctx); printf("\n");
                flint_printf("lambda2: "); ca_vec_print(lambda2, ctx); printf("\n");
                flint_abort();
            }
        }
        else
        {
            flint_printf("FAIL (unexpected failure to compute Jordan blocks in rational case)\n");
            flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("lambda1: "); ca_vec_print(lambda1, ctx); printf("\n");
            flint_abort();
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(P, ctx);
        ca_mat_clear(Q, ctx);
        ca_mat_clear(B, ctx);
        ca_vec_clear(lambda1, ctx);
        ca_vec_clear(lambda2, ctx);

        flint_free(block_lambda1);
        flint_free(block_lambda2);
        flint_free(block_size1);
        flint_free(block_size2);

        ca_ctx_clear(ctx);
    }

    /* Test J(A) = J(P * A * P^-1) for random A, P */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, P, Q, B;
        slong n;
        int success1, success2;
        slong num_blocks1, num_blocks2;
        ca_vec_t lambda1, lambda2;
        slong *block_lambda1, *block_size1, *block_lambda2, *block_size2;

        n = n_randint(state, 6);

        ca_ctx_init(ctx);
        ca_mat_init(A, n, n, ctx);
        ca_mat_init(P, n, n, ctx);
        ca_mat_init(Q, n, n, ctx);
        ca_mat_init(B, n, n, ctx);
        ca_vec_init(lambda1, 0, ctx);
        ca_vec_init(lambda2, 0, ctx);

        block_lambda1 = flint_malloc(sizeof(slong) * n);
        block_lambda2 = flint_malloc(sizeof(slong) * n);
        block_size1 = flint_malloc(sizeof(slong) * n);
        block_size2 = flint_malloc(sizeof(slong) * n);

        if (n <= 3 && n_randint(state, 2) == 0)
            ca_mat_randtest(A, state, 1, 5, ctx);
        else
            ca_mat_randtest_rational(A, state, 5, ctx);

        do {
            ca_mat_randtest_rational(P, state, 1, ctx);
        } while (ca_mat_inv(Q, P, ctx) != T_TRUE);

        ca_mat_mul(B, P, A, ctx);
        ca_mat_mul(B, B, Q, ctx);

        success1 = ca_mat_jordan_blocks(lambda1, &num_blocks1, block_lambda1, block_size1, A, ctx);
        success2 = success1 && ca_mat_jordan_blocks(lambda2, &num_blocks2, block_lambda2, block_size2, B, ctx);

        if (success1 && success2)
        {
            if (!check_jordan_forms(lambda1, num_blocks1, block_lambda1, block_size1,
                                    lambda2, num_blocks2, block_lambda2, block_size2, ctx))
            {
                flint_printf("FAIL (different jordan forms)\n");
                flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("B: "); ca_mat_print(B, ctx); flint_printf("\n");
                flint_printf("lambda1: "); ca_vec_print(lambda1, ctx); printf("\n");
                flint_printf("lambda2: "); ca_vec_print(lambda2, ctx); printf("\n");
                flint_abort();
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(P, ctx);
        ca_mat_clear(Q, ctx);
        ca_mat_clear(B, ctx);
        ca_vec_clear(lambda1, ctx);
        ca_vec_clear(lambda2, ctx);

        flint_free(block_lambda1);
        flint_free(block_lambda2);
        flint_free(block_size1);
        flint_free(block_size2);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
