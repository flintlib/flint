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

int test_conversion(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    slong M = 20;
    slong N = 20;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    gr_coo_mat_t coo_mat, coo_mat2;
    gr_csr_mat_t csr_mat, csr_mat2;
    gr_lil_mat_t lil_mat, lil_mat2;
    gr_mat_t dmat, dmat2;

    gr_mat_init(dmat, M, N, ctx);
    gr_mat_init(dmat2, M, N, ctx);
    gr_coo_mat_init(coo_mat, M, N, ctx);
    gr_coo_mat_init(coo_mat2, M, N, ctx);
    gr_csr_mat_init(csr_mat, M, N, ctx);
    gr_csr_mat_init(csr_mat2, M, N, ctx);
    gr_lil_mat_init(lil_mat, M, N, ctx);
    gr_lil_mat_init(lil_mat2, M, N, ctx);

    //flint_printf("Testing copy\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_coo_mat_randtest(coo_mat, 40, 0, T_TRUE, state, ctx);
        status |= gr_coo_mat_set(coo_mat2, coo_mat, ctx);
        if (gr_coo_mat_equal(coo_mat, coo_mat2, ctx) == T_FALSE)
        {
            status |= gr_coo_mat_print_nz(coo_mat, ctx); flint_printf("\n");
            status |= gr_coo_mat_print_nz(coo_mat2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    //flint_printf("Testing from/to dense mat\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_mat_randtest(dmat, state, ctx);
        status |= gr_coo_mat_set_mat(coo_mat, dmat, ctx);
        status |= gr_mat_set_coo_mat(dmat2, coo_mat, ctx);
        if (T_FALSE == gr_mat_equal(dmat, dmat2, ctx))
        {
            status |= gr_mat_print(dmat, ctx); flint_printf("\n");
            status |= gr_mat_print(dmat2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    // flint_printf("Testing from/to sparse mat\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_coo_mat_randtest(coo_mat, 10, 0, T_TRUE, state, ctx);
        status |= gr_mat_set_coo_mat(dmat, coo_mat, ctx);
        status |= gr_coo_mat_set_mat(coo_mat2, dmat, ctx);
        if (T_FALSE == gr_coo_mat_equal(coo_mat, coo_mat2, ctx))
        {
            status |= gr_coo_mat_print_nz(coo_mat, ctx); flint_printf("\n");
            status |= gr_coo_mat_print_nz(coo_mat2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    // flint_printf("Testing coo -> csr -> csr -> lil -> mat -> csr -> coo \n");
    for (i = 0; i < 2 * n_tests; i++)
    {
        status |= gr_coo_mat_randtest(coo_mat, 10, 0, (i % 2) ? T_FALSE : T_TRUE, state, ctx);
        //flint_printf("\n\ncoo_mat = "); status |= gr_coo_mat_print_nz(coo_mat, ctx);
        status |= gr_csr_mat_set_coo_mat(csr_mat, coo_mat, ctx);
        //flint_printf("\n\ncsr_mat = "); status |= gr_csr_mat_print_nz(csr_mat, ctx); flint_printf("\nnnz = %d\n", csr_mat->nnz);
        status |= gr_csr_mat_set(csr_mat2, csr_mat, ctx);
        //flint_printf("\n\ncsr_mat = "); status |= gr_csr_mat_print_nz(csr_mat2, ctx); flint_printf("\nnnz = %d\n", csr_mat2->nnz);
        status |= gr_lil_mat_set_csr_mat(lil_mat, csr_mat2, ctx);
        //flint_printf("\n\nlil_mat = "); status |= gr_lil_mat_print_nz(lil_mat, ctx);
        status |= gr_mat_set_lil_mat(dmat, lil_mat, ctx);
        //flint_printf("\n\nmat = "); status |= gr_mat_print(dmat, ctx);
        status |= gr_csr_mat_set_mat(csr_mat2, dmat, ctx);
        //flint_printf("\n\ncsr_mat = "); status |= gr_csr_mat_print_nz(csr_mat2, ctx);
        status |= gr_coo_mat_set_csr_mat(coo_mat2, csr_mat2, ctx);
        //flint_printf("\n\ncoo_mat = "); status |= gr_coo_mat_print_nz(coo_mat2, ctx); flint_printf("\nnnz = %d\n", coo_mat2->nnz);
        if 
        (
            T_FALSE == gr_csr_mat_equal(csr_mat, csr_mat2, ctx) ||
            T_FALSE == gr_coo_mat_equal(coo_mat, coo_mat2, ctx)
        )
        {
            status |= gr_coo_mat_print_nz(coo_mat, ctx); flint_printf("\n");
            status |= gr_csr_mat_print_nz(csr_mat2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    // flint_printf("Testing coo -> lil -> lil -> csr -> mat -> lil -> coo\n");
    for (i = 0; i < 2 * n_tests; i++)
    {
        status |= gr_coo_mat_randtest(coo_mat, 10, 0, (i % 2) ? T_FALSE : T_TRUE, state, ctx);
        // flint_printf("\n\ncoo_mat = "); status |= gr_coo_mat_print_nz(coo_mat, ctx); flint_printf("\nnnz = %d\n", coo_mat->nnz);
        status |= gr_lil_mat_set_coo_mat(lil_mat, coo_mat, ctx);
        // flint_printf("\n\nlil_mat = "); status |= gr_lil_mat_print_nz(lil_mat, ctx); flint_printf("\nnnz = %d\n", lil_mat->nnz);
        status |= gr_lil_mat_set(lil_mat2, lil_mat, ctx);
        // flint_printf("\n\nlil_mat = "); status |= gr_lil_mat_print_nz(lil_mat2, ctx); flint_printf("\nnnz = %d\n", lil_mat2->nnz);
        status |= gr_csr_mat_set_lil_mat(csr_mat, lil_mat2, ctx);
        // flint_printf("\n\ncsr_mat = "); status |= gr_csr_mat_print_nz(csr_mat, ctx); flint_printf("\nnnz = %d\n", csr_mat->nnz);
        status |= gr_mat_set_csr_mat(dmat, csr_mat, ctx);
        // flint_printf("\n\nmat = "); status |= gr_mat_print(dmat, ctx);
        status |= gr_lil_mat_set_mat(lil_mat2, dmat, ctx);
        // flint_printf("\n\nlil_mat = "); status |= gr_lil_mat_print_nz(lil_mat2, ctx); flint_printf("\nnnz = %d\n", lil_mat2->nnz);
        status |= gr_coo_mat_set_lil_mat(coo_mat2, lil_mat2, ctx);
        // flint_printf("\n\ncoo_mat = "); status |= gr_coo_mat_print_nz(coo_mat2, ctx); flint_printf("\nnnz = %d\n", coo_mat2->nnz);
        if
        (
            T_FALSE == gr_lil_mat_equal(lil_mat, lil_mat2, ctx) ||
            T_FALSE == gr_coo_mat_equal(coo_mat, coo_mat2, ctx)
        )
        {
            status |= gr_coo_mat_print_nz(coo_mat, ctx); flint_printf("\n");
            status |= gr_lil_mat_print_nz(lil_mat2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    gr_csr_mat_clear(csr_mat, ctx);
    gr_csr_mat_clear(csr_mat2, ctx);
    gr_lil_mat_clear(lil_mat, ctx);
    gr_lil_mat_clear(lil_mat2, ctx);
    gr_coo_mat_clear(coo_mat, ctx);
    gr_coo_mat_clear(coo_mat2, ctx);
    gr_mat_clear(dmat, ctx);
    gr_mat_clear(dmat2, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_mat_conversion, state)
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
        CHECK_TEST(test_conversion(state, ctx), "Conversion between various sparse representations (and dense)");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
