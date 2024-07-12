/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_sparse_mat.h"

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }

int test_randtest(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    slong M = 16;
    slong N = 64;
    slong n_tests = 10;
    int status = GR_SUCCESS;
    gr_coo_mat_t mat;
    gr_coo_mat_init(mat, M, N, ctx);

    //flint_printf("Testing w/o replacement\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_coo_mat_randtest(mat, 128, 0, T_TRUE, state, ctx);
        if (gr_coo_mat_is_canonical(mat, ctx) == T_FALSE || mat->nnz != 128)
            return GR_TEST_FAIL;
    }

    //flint_printf("Testing w/ replacement\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_coo_mat_randtest(mat, 32, 1, T_TRUE, state, ctx);
        if (gr_coo_mat_is_canonical(mat, ctx) == T_FALSE || mat->nnz > 32 || mat->nnz < 24)
            return GR_TEST_FAIL;
    }

    //flint_printf("Testing w/ prob\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_coo_mat_randtest_prob(mat, 0.125, state, ctx);
        if (gr_coo_mat_is_canonical(mat, ctx) == T_FALSE || mat->nnz > 192 || mat->nnz < 64)
        {
            status |= gr_coo_mat_print_nz(mat, ctx); flint_printf("%ld\n", mat->nnz);
            return GR_TEST_FAIL;
        }
            
    }
    gr_coo_mat_clear(mat, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_mat_randtest, state)
{   
    int i;
    gr_ctx_t ctx;
    for (i = 0; i < 16; ++i)
    {
        while(1)
        {
            //gr_ctx_init_nmod(ctx, 2147483647);
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        //gr_ctx_print(ctx); flint_printf("\n");

        CHECK_TEST(test_randtest(state, ctx), "Test random sparse matrix generation");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
