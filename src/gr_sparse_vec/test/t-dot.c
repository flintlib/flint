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

int test_dot(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    slong N = 30;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    gr_sparse_vec_t vec, vec2;
    gr_ptr dvec, dvec2;
    gr_ptr dot, dot2;
    GR_TMP_INIT2(dot, dot2, ctx);

    GR_TMP_INIT_VEC(dvec, N, ctx);
    GR_TMP_INIT_VEC(dvec2, N, ctx);
    gr_sparse_vec_init(vec, N, ctx);
    gr_sparse_vec_init(vec2, N, ctx);

    //flint_printf("Testing copy\n");
    for (i = 0; i < 2*n_tests; i++)
    {
        // Give a random initial value to dot
        status |= gr_randtest(dot, state, ctx);
        status |= gr_set(dot2, dot, ctx);

        // Get two random sparse vectors
        status |= gr_sparse_vec_randtest(vec, 10, 0, state, ctx);
        status |= gr_sparse_vec_randtest(vec2, 10, 0, state, ctx);
        status |= gr_vec_set_sparse_vec(dvec, vec, ctx);
        status |= gr_vec_set_sparse_vec(dvec2, vec2, ctx);

        // Compute sparse dot product and check against dense one
        status |= gr_sparse_vec_dot(dot, dot, i % 2, vec, vec2, ctx);
        status |= _gr_vec_dot(dot2, dot2, i % 2, dvec, dvec2, N, ctx);
        if (gr_equal(dot, dot2, ctx) == T_FALSE)
        {
            gr_ctx_println(ctx);
            flint_printf("vec = "); gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
            flint_printf("vec2 = "); gr_sparse_vec_print_nz(vec2, ctx); flint_printf("\n");
            flint_printf("dvec = "); _gr_vec_print(dvec, N, ctx); flint_printf("\n");
            flint_printf("dvec2 = "); _gr_vec_print(dvec2, N, ctx); flint_printf("\n");
            flint_printf("dot = "); gr_println(dot, ctx);
            flint_printf("dot2 = "); gr_println(dot2, ctx);
            return GR_TEST_FAIL;
        }
    }

    //flint_printf("Testing from/to dense vec\n");
    for (i = 0; i < 2*n_tests; i++)
    {
        // Give a random initial value to dot
        status |= gr_randtest(dot, state, ctx);
        status |= gr_set(dot2, dot, ctx);

        // Get random sparse and dense vectors
        status |= gr_sparse_vec_randtest(vec, 10, 0, state, ctx);
        status |= _gr_vec_randtest(dvec2, state, N, ctx);
        status |= gr_vec_set_sparse_vec(dvec, vec, ctx);

        // Compute sparse-dense dot product and check against dense-dense one
        status |= gr_sparse_vec_dot_vec(dot, dot, i % 2, vec, dvec2, ctx);
        status |= _gr_vec_dot(dot2, dot2, i % 2, dvec, dvec2, N, ctx);

        if (gr_equal(dot, dot2, ctx) == T_FALSE)
        {
            gr_ctx_println(ctx);
            flint_printf("vec = "); gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
            flint_printf("vec2 = "); _gr_vec_print(dvec2, N, ctx); flint_printf("\n");
            flint_printf("dot = "); gr_println(dot, ctx);
            flint_printf("dot2 = "); gr_println(dot2, ctx);
            return GR_TEST_FAIL;
        }
    }

    gr_sparse_vec_clear(vec, ctx);
    gr_sparse_vec_clear(vec2, ctx);
    GR_TMP_CLEAR_VEC(dvec, N, ctx);
    GR_TMP_CLEAR_VEC(dvec2, N, ctx);
    GR_TMP_CLEAR2(dot, dot2, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_vec_dot, state)
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
        CHECK_TEST(test_dot(state, ctx), "Dot product");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
