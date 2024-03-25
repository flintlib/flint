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

int test_mul(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j;
    slong M = 5;
    slong N = 15;
    slong O = 2;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    truth_t eq;
    gr_csr_mat_t csr_mat;
    gr_lil_mat_t lil_mat;
    gr_coo_mat_t coo_mat;
    gr_ptr u, v, v2;
    gr_mat_t dmat, U, V, V2;
    gr_mat_t UT, VT, VT2;

    gr_coo_mat_init(coo_mat, M, N, ctx);
    gr_csr_mat_init(csr_mat, M, N, ctx);
    gr_lil_mat_init(lil_mat, M, N, ctx);
    gr_mat_init(dmat, M, N, ctx);

    GR_TMP_INIT_VEC(u, N, ctx);
    gr_mat_init(UT, O, N, ctx);
    gr_mat_init(U, N, O, ctx);

    GR_TMP_INIT_VEC(v, M, ctx);
    gr_mat_init(VT, O, M, ctx);
    gr_mat_init(V, M, O, ctx);

    GR_TMP_INIT_VEC(v2, M, ctx);
    gr_mat_init(VT2, O, M, ctx);
    gr_mat_init(V2, M, O, ctx);

    //flint_printf("Testing sparse matrix-vector multiplication\n");
    for (i = 0; i < 2*n_tests; i++)
    {
        j = i % 2; // CSR or LIL mat

        // Get random sparse matrix
        status |= gr_coo_mat_randtest(coo_mat, 20, 0, T_TRUE, state, ctx);
        status |= gr_mat_set_coo_mat(dmat, coo_mat, ctx);
        if (j == 0)
            status |= gr_csr_mat_set_coo_mat(csr_mat, coo_mat, ctx);
        else
            status |= gr_lil_mat_set_coo_mat(lil_mat, coo_mat, ctx);

        // Get a random dense vector
        status |= _gr_vec_randtest(u, state, N, ctx);

        // Compute matrix multiply using sparse and dense reprs
        if (j == 0)
            status |= gr_csr_mat_mul_vec(v, csr_mat, u, ctx);
        else
            status |= gr_lil_mat_mul_vec(v, lil_mat, u, ctx);
        status |= gr_mat_mul_vec(v2, dmat, u, ctx);
        
        eq = _gr_vec_equal(v, v2, M, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            gr_ctx_println(ctx);
            if (j == 0)
            {
                status |= flint_printf("mat = "); status |= gr_csr_mat_print_nz(csr_mat, ctx); flint_printf("\n");
            }
            else
            {
                status |= flint_printf("mat = "); status |= gr_lil_mat_print_nz(lil_mat, ctx); flint_printf("\n");
            }
            flint_printf("u = "); status |= _gr_vec_print(u, N, ctx); flint_printf("\n");
            flint_printf("v = "); status |= _gr_vec_print(v, M, ctx); flint_printf("\n");
            flint_printf("v2 = "); status |= _gr_vec_print(v2, M, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }

        // Get a random dense column matrix
        status |= gr_mat_randtest(UT, state, ctx);

        // Compute matrix multiply using sparse and dense reprs
        if (j == 0)
            status |= gr_csr_mat_mul_mat_transpose(VT, csr_mat, UT, ctx);
        else
            status |= gr_lil_mat_mul_mat_transpose(VT, lil_mat, UT, ctx);
        status |= gr_mat_transpose(U, UT, ctx);
        status |= gr_mat_mul(V2, dmat, U, ctx);
        status |= gr_mat_transpose(VT2, V2, ctx);
                
        eq = gr_mat_equal(VT, VT2, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            gr_ctx_println(ctx);
            flint_printf("j = %d, eq = %d, status = %d\n", j, eq, status);
            if (j == 0)
            {
                status |= flint_printf("mat = "); status |= gr_csr_mat_print_nz(csr_mat, ctx); flint_printf("\n");
            }
            else
            {
                status |= flint_printf("mat = "); status |= gr_lil_mat_print_nz(lil_mat, ctx); flint_printf("\n");
            }
            flint_printf("UT = "); status |= gr_mat_print(UT, ctx); flint_printf("\n");
            flint_printf("VT = "); status |= gr_mat_print(VT, ctx); flint_printf("\n");
            flint_printf("VT2 = "); status |= gr_mat_print(VT2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }

        // Get a random dense row matrix
        status |= gr_mat_randtest(U, state, ctx);

        // Compute matrix multiply using sparse and dense reprs
        if (j == 0)
            status |= gr_csr_mat_mul_mat(V, csr_mat, U, ctx);
        else
            status |= gr_lil_mat_mul_mat(V, lil_mat, U, ctx);

        //flint_printf("\nBefore matmul U = "); gr_mat_print(U, ctx); flint_printf("\n");
        status |= gr_mat_mul(V2, dmat, U, ctx);
        //flint_printf("\nAfter matmul U = "); gr_mat_print(U, ctx); flint_printf("\n");
                
        eq = gr_mat_equal(V, V2, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            gr_ctx_println(ctx);
            flint_printf("j = %d, eq = %d, status = %d\n", j, eq, status);
            if (j == 0)
            {
                status |= flint_printf("mat = "); status |= gr_csr_mat_print_nz(csr_mat, ctx); flint_printf("\n");
            }
            else
            {
                status |= flint_printf("mat = "); status |= gr_lil_mat_print_nz(lil_mat, ctx); flint_printf("\n");
            }
            flint_printf("U = "); status |= gr_mat_print(U, ctx); flint_printf("\n");
            flint_printf("V = "); status |= gr_mat_print(V, ctx); flint_printf("\n");
            flint_printf("V2 = "); status |= gr_mat_print(V2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    gr_mat_clear(U, ctx);
    gr_mat_clear(UT, ctx);
    GR_TMP_CLEAR_VEC(u, N, ctx);

    gr_mat_clear(V, ctx);
    gr_mat_clear(VT, ctx);
    GR_TMP_CLEAR_VEC(v, M, ctx);

    gr_mat_clear(VT2, ctx);
    gr_mat_clear(V2, ctx);
    GR_TMP_CLEAR_VEC(v2, M, ctx);

    gr_mat_clear(dmat, ctx);
    gr_csr_mat_clear(csr_mat, ctx);
    gr_lil_mat_clear(lil_mat, ctx);
    gr_coo_mat_clear(coo_mat, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_mat_mul, state)
{   
    slong i;
    gr_ctx_t ctx;
    for (i = 0; i < 16; ++i)
    {
        while (1)
        {
            //gr_ctx_init_fmpz(ctx);
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        //gr_ctx_println(ctx);
        CHECK_TEST(test_mul(state, ctx), "Sparse matrix-vector and sparse matrix-matrix products");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
