/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_sparse_vec.h"

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }

int test_randtest(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    slong N = 1024;
    slong n_tests = 10;
    int status = GR_SUCCESS;
    gr_sparse_vec_t vec;
    gr_sparse_vec_init(vec, N, ctx);

    //flint_printf("Testing w/o replacement\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_sparse_vec_randtest(vec, 128, 0, state, ctx);
        if (!gr_sparse_vec_is_valid(vec, ctx) || vec->nnz != 128)
            return GR_TEST_FAIL;
    }

    //flint_printf("Testing w/ replacement\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_sparse_vec_randtest(vec, 32, 1, state, ctx);
        if (!gr_sparse_vec_is_valid(vec, ctx) || vec->nnz > 32 || vec->nnz < 24)
            return GR_TEST_FAIL;
    }

    //flint_printf("Testing w/ prob\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_sparse_vec_randtest_prob(vec, 0.125, state, ctx);
        if (!gr_sparse_vec_is_valid(vec, ctx) || vec->nnz > 192 || vec->nnz < 64)
        {
            gr_sparse_vec_print_nz(vec, ctx); flint_printf("%ld\n", vec->nnz);
            return GR_TEST_FAIL;
        }
            
    }
    gr_sparse_vec_clear(vec, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_vec_randtest, state)
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

        CHECK_TEST(test_randtest(state, ctx), "Test random sparse vector generation");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
