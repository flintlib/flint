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

int test_sum_prod(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    slong N = 30;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    gr_sparse_vec_t vec;
    gr_ptr dvec;
    gr_ptr res, res2;
    GR_TMP_INIT2(res, res2, ctx);

    GR_TMP_INIT_VEC(dvec, N, ctx);
    gr_sparse_vec_init(vec, N, ctx);

    //flint_printf("Testing sum\n");
    for (i = 0; i < 2*n_tests; i++)
    {
        // Get random sparse vector
        status |= gr_sparse_vec_randtest(vec, 10, 0, state, ctx);
        status |= gr_vec_set_sparse_vec(dvec, vec, ctx);

        // Compute sparse sum and check against dense one
        status |= gr_sparse_vec_sum(res, vec, ctx);
        status |= _gr_vec_sum(res2, dvec, N, ctx);
        if (gr_equal(res, res2, ctx) == T_FALSE)
        {
            gr_ctx_println(ctx);
            flint_printf("vec = "); gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
            flint_printf("res = "); gr_println(res, ctx);
            flint_printf("res2 = "); gr_println(res2, ctx);
            return GR_TEST_FAIL;
        }
    }

    //flint_printf("Testing prod\n");
    for (i = 0; i < 2*n_tests; i++)
    {
        // Get random sparse vector
        status |= gr_sparse_vec_randtest(vec, 10, 0, state, ctx);
        status |= gr_vec_set_sparse_vec(dvec, vec, ctx);

        // Compute sparse sum and check against dense one
        status |= gr_sparse_vec_nz_product(res, vec, ctx);
        status |= _gr_vec_product(res2, vec->nzs, vec->nnz, ctx);
        if (gr_equal(res, res2, ctx) == T_FALSE)
        {
            gr_ctx_println(ctx);
            flint_printf("vec = "); gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
            flint_printf("res = "); gr_println(res, ctx);
            flint_printf("res2 = "); gr_println(res2, ctx);
            return GR_TEST_FAIL;
        }
    }

    gr_sparse_vec_clear(vec, ctx);
    GR_TMP_CLEAR_VEC(dvec, N, ctx);
    GR_TMP_CLEAR2(res, res2, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_vec_sum_prod, state)
{   
    slong i;
    gr_ctx_t ctx;
    for (i = 0; i < 16; ++i)
    {
        while (1)
        {
            gr_ctx_init_fmpz(ctx);
            //gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        //gr_ctx_println(ctx);
        CHECK_TEST(test_sum_prod(state, ctx), "Sum and (nonzero) product");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
