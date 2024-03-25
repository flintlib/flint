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
#include "fmpz.h"
#include "fmpq.h"

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }

int test_neg(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    slong M = 20;
    slong N = 20;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    gr_csr_mat_t csr_mat, csr_mat2;
    gr_lil_mat_t lil_mat, lil_mat2;
    gr_coo_mat_t coo_mat, coo_mat2;
    gr_mat_t dmat, dmat2;

    gr_mat_init(dmat, M, N, ctx);
    gr_mat_init(dmat2, M, N, ctx);
    gr_coo_mat_init(coo_mat, M, N, ctx);
    gr_coo_mat_init(coo_mat2, M, N, ctx);
    gr_csr_mat_init(csr_mat, M, N, ctx);
    gr_csr_mat_init(csr_mat2, M, N, ctx);
    gr_lil_mat_init(lil_mat, M, N, ctx);
    gr_lil_mat_init(lil_mat2, M, N, ctx);

    //flint_printf("Testing neg coo\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_coo_mat_randtest(coo_mat, 80, 0, T_TRUE, state, ctx);
        status |= gr_coo_mat_neg(coo_mat2, coo_mat, ctx);
        status |= gr_mat_set_coo_mat(dmat2, coo_mat, ctx);
        status |= gr_mat_neg(dmat, dmat2, ctx);
        status |= gr_mat_set_coo_mat(dmat2, coo_mat2, ctx);
        if (gr_mat_equal(dmat, dmat2, ctx) == T_FALSE)
        {
            status |= gr_coo_mat_print_nz(coo_mat, ctx); flint_printf("\n");
            status |= gr_coo_mat_print_nz(coo_mat2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    //flint_printf("Testing neg csr\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_coo_mat_randtest(coo_mat, 80, 0, T_TRUE, state, ctx);
        status |= gr_csr_mat_set_coo_mat(csr_mat, coo_mat, ctx);
        status |= gr_csr_mat_neg(csr_mat2, csr_mat, ctx);
        status |= gr_mat_set_csr_mat(dmat2, csr_mat, ctx);
        status |= gr_mat_neg(dmat, dmat2, ctx);
        status |= gr_mat_set_csr_mat(dmat2, csr_mat2, ctx);
        if (gr_mat_equal(dmat, dmat2, ctx) == T_FALSE)
        {
            status |= gr_csr_mat_print_nz(csr_mat, ctx); flint_printf("\n");
            status |= gr_csr_mat_print_nz(csr_mat2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    //flint_printf("Testing neg lil\n");
    for (i = 0; i < n_tests; i++)
    {
        status |= gr_coo_mat_randtest(coo_mat, 80, 0, T_TRUE, state, ctx);
        status |= gr_lil_mat_set_coo_mat(lil_mat, coo_mat, ctx);
        status |= gr_lil_mat_neg(lil_mat2, lil_mat, ctx);
        status |= gr_mat_set_lil_mat(dmat2, lil_mat, ctx);
        status |= gr_mat_neg(dmat, dmat2, ctx);
        status |= gr_mat_set_lil_mat(dmat2, lil_mat2, ctx);
        if (gr_mat_equal(dmat, dmat2, ctx) == T_FALSE)
        {
            status |= gr_lil_mat_print_nz(lil_mat, ctx); flint_printf("\n");
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

int test_add_sub(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j;
    slong M = 20;
    slong N = 20;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    truth_t eq;
    gr_csr_mat_t csr_mat, csr_mat2;
    gr_lil_mat_t lil_mat, lil_mat2;
    gr_coo_mat_t coo_mat, coo_mat2;
    gr_mat_t dmat, dmat2;

    gr_mat_init(dmat, M, N, ctx);
    gr_mat_init(dmat2, M, N, ctx);
    gr_coo_mat_init(coo_mat, M, N, ctx);
    gr_coo_mat_init(coo_mat2, M, N, ctx);
    gr_csr_mat_init(csr_mat, M, N, ctx);
    gr_csr_mat_init(csr_mat2, M, N, ctx);
    gr_lil_mat_init(lil_mat, M, N, ctx);
    gr_lil_mat_init(lil_mat2, M, N, ctx);

    //flint_printf("Testing add sub\n");
    for (i = 0; i < 2 * n_tests; i++)
    {
        status = GR_SUCCESS;
        //flint_printf("%d\n", i);
        j = i % 2; // Add or subtract

        status |= gr_coo_mat_randtest(coo_mat, 80, 0, T_TRUE, state, ctx);
        status |= gr_coo_mat_randtest(coo_mat2, 80, 0, T_TRUE, state, ctx);
        status |= gr_lil_mat_set_coo_mat(lil_mat, coo_mat, ctx);
        status |= gr_lil_mat_set_coo_mat(lil_mat2, coo_mat2, ctx);
        status |= gr_mat_set_lil_mat(dmat, lil_mat, ctx);
        status |= gr_mat_set_lil_mat(dmat2, lil_mat2, ctx);
        if (j == 0)
        {
            status |= gr_lil_mat_add(lil_mat, lil_mat, lil_mat2, ctx);
            status |= gr_mat_add(dmat, dmat, dmat2, ctx);
        }
        else
        {
            status |= gr_lil_mat_sub(lil_mat, lil_mat, lil_mat2, ctx);
            status |= gr_mat_sub(dmat, dmat, dmat2, ctx);
        }

        if (status == GR_UNABLE)
            continue;

        status |= gr_mat_set_lil_mat(dmat2, lil_mat, ctx);
        eq = gr_mat_equal(dmat, dmat2, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            flint_printf(
                "\ni = %d, j = %d, equal = %d, status = %d\n",
                i, j, eq, status
            );
            gr_ctx_println(ctx);
            flint_printf("mat = "); gr_mat_print(dmat, ctx); flint_printf("\n");
            flint_printf("mat2 = "); gr_mat_print(dmat2, ctx); flint_printf("\n");
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

int test_accum_mul_scalar(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j;
    slong M = 20;
    slong N = 20;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    truth_t eq;
    gr_lil_mat_t lil_mat, lil_mat2;
    gr_coo_mat_t coo_mat, coo_mat2;
    gr_mat_t dmat, dmat2;
    gr_ptr gr_c;

    GR_TMP_INIT(gr_c, ctx);
    gr_mat_init(dmat, M, N, ctx);
    gr_mat_init(dmat2, M, N, ctx);
    gr_coo_mat_init(coo_mat, M, N, ctx);
    gr_coo_mat_init(coo_mat2, M, N, ctx);
    gr_lil_mat_init(lil_mat, M, N, ctx);
    gr_lil_mat_init(lil_mat2, M, N, ctx);

    //flint_printf("Testing addmul submul\n");
    for (i = 0; i < 2 * n_tests; i++)
    {
        status = GR_SUCCESS;
        //flint_printf("%d\n", i);
        j = i % 2; // Add or subtract

        status |= gr_coo_mat_randtest(coo_mat, 80, 0, T_TRUE, state, ctx);
        status |= gr_coo_mat_randtest(coo_mat2, 80, 0, T_TRUE, state, ctx);
        status |= gr_lil_mat_set_coo_mat(lil_mat, coo_mat, ctx);
        status |= gr_lil_mat_set_coo_mat(lil_mat2, coo_mat2, ctx);
        status |= gr_mat_set_lil_mat(dmat, lil_mat, ctx);
        status |= gr_mat_set_lil_mat(dmat2, lil_mat2, ctx);

        status |= gr_randtest_not_zero(gr_c, state, ctx);
        if (j == 0)
        {
            status |= gr_lil_mat_addmul_scalar(lil_mat, lil_mat2, gr_c, ctx);
            status |= gr_mat_addmul_scalar(dmat, dmat2, gr_c, ctx);
        }
        else
        {
            status |= gr_lil_mat_submul_scalar(lil_mat, lil_mat2, gr_c, ctx);
            status |= gr_mat_submul_scalar(dmat, dmat2, gr_c, ctx);
        }
        status |= gr_mat_set_lil_mat(dmat2, lil_mat, ctx);
        eq = gr_mat_equal(dmat, dmat2, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            flint_printf(
                "\ni = %d, j = %d, equal = %d, status = %d\n",
                i, j, eq, status
            );
            gr_ctx_println(ctx);
            flint_printf("dmat = "); gr_mat_print(dmat, ctx); flint_printf("\n");
            flint_printf("dmat2 = "); gr_mat_print(dmat2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    gr_lil_mat_clear(lil_mat, ctx);
    gr_lil_mat_clear(lil_mat2, ctx);
    gr_coo_mat_clear(coo_mat, ctx);
    gr_coo_mat_clear(coo_mat2, ctx);
    gr_mat_clear(dmat, ctx);
    gr_mat_clear(dmat2, ctx);
    GR_TMP_CLEAR(gr_c, ctx);
    return status;
}

#define TEST_OP_SCALAR(STATUS, MAT_TYPE, OP, SCALAR_TYPE, MAT, MAT2, DMAT, DMAT2, C, CTX) \
        STATUS |= gr_##MAT_TYPE##_mat_##OP##_##SCALAR_TYPE(MAT2, MAT, C, CTX); \
        STATUS |= gr_mat_##OP##_##SCALAR_TYPE(DMAT2, DMAT, C, CTX); \


#define TEST_MUL_SCALAR(STATUS, K, MAT_TYPE, SCALAR_TYPE, MAT, MAT2, DMAT, DMAT2, C, CTX) { \
    if (K == 1) \
    { TEST_OP_SCALAR(STATUS, MAT_TYPE, div, SCALAR_TYPE, MAT, MAT2, DMAT, DMAT2, C, CTX) } \
    else \
    { \
        TEST_OP_SCALAR(STATUS, MAT_TYPE, mul, SCALAR_TYPE, MAT, MAT2, DMAT, DMAT2, C, CTX) \
        if (K == 2) \
        { TEST_OP_SCALAR(STATUS, MAT_TYPE, divexact, SCALAR_TYPE, MAT2, MAT2, DMAT2, DMAT2, C, CTX) } \
    } \
}

int test_mul_div_scalar(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j, k, l;
    slong M = 20;
    slong N = 20;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    gr_csr_mat_t csr_mat, csr_mat2;
    gr_lil_mat_t lil_mat, lil_mat2;
    gr_coo_mat_t coo_mat, coo_mat2;
    gr_mat_t dmat, dmat2;
    slong c;
    ulong uc;
    gr_ptr gr_c;
    fmpz_t fmpz_c;
    fmpq_t fmpq_c;
    truth_t eq;

    fmpz_init(fmpz_c);
    fmpq_init(fmpq_c);

    GR_TMP_INIT(gr_c, ctx);
    gr_mat_init(dmat, M, N, ctx);
    gr_mat_init(dmat2, M, N, ctx);
    gr_coo_mat_init(coo_mat, M, N, ctx);
    gr_coo_mat_init(coo_mat2, M, N, ctx);
    gr_csr_mat_init(csr_mat, M, N, ctx);
    gr_csr_mat_init(csr_mat2, M, N, ctx);
    gr_lil_mat_init(lil_mat, M, N, ctx);
    gr_lil_mat_init(lil_mat2, M, N, ctx);

    //flint_printf("Testing mul div scalar\n");
    for (i = 0; i < 54 * n_tests; i++)
    {
        j = i % 6; // Which type of scalar
        k = (i / 6) % 3; // Mul, div, or mul + divexact
        l = (i / 18) % 3; // CSR, LIL, or COO mat
        //flint_printf("\nTesting (%d, %d, %d)\n", l, k, j);
        if ((j == 4 || k == 1) && gr_ctx_is_field(ctx) != T_TRUE)
            continue;
        if (k == 2 && (gr_ctx_is_exact(ctx) != T_TRUE || gr_ctx_is_integral_domain(ctx) != T_TRUE))
            continue;
        status |= gr_coo_mat_randtest(coo_mat, 80, 0, T_TRUE, state, ctx);
        status |= gr_mat_set_coo_mat(dmat, coo_mat, ctx);
        if (l == 0)
            status |= gr_csr_mat_set_coo_mat(csr_mat, coo_mat, ctx);
        else if (l == 1)
            status |= gr_lil_mat_set_coo_mat(lil_mat, coo_mat, ctx);

        //flint_printf("\nmat = "); status |= gr_coo_mat_print_nz(coo_mat, ctx); flint_printf("\n");
        //flint_printf("\nmat = "); status |= gr_mat_print(dmat, ctx); flint_printf("\n");
        switch(j)
        {
        case 0: 
            //flint_printf("Testing scalar\n");
            status |= gr_randtest_not_zero(gr_c, state, ctx);
            if (l == 0)
                TEST_MUL_SCALAR(status, k, csr, scalar, csr_mat, csr_mat2, dmat, dmat2, gr_c, ctx)
            else if (l == 1)
                TEST_MUL_SCALAR(status, k, lil, scalar, lil_mat, lil_mat2, dmat, dmat2, gr_c, ctx)
            else
                TEST_MUL_SCALAR(status, k, coo, scalar, coo_mat, coo_mat2, dmat, dmat2, gr_c, ctx)
            break;
        case 1:
            c = n_randint(state, 0);
            //flint_printf("Testing scalar_si = %ld\n", c);
            if (l == 0)
                TEST_MUL_SCALAR(status, k, csr, scalar_si, csr_mat, csr_mat2, dmat, dmat2, c, ctx)
            else if (l == 1)
                TEST_MUL_SCALAR(status, k, lil, scalar_si, lil_mat, lil_mat2, dmat, dmat2, c, ctx)
            else
                TEST_MUL_SCALAR(status, k, coo, scalar_si, coo_mat, coo_mat2, dmat, dmat2, c, ctx)
            break;
        case 2:
            //flint_printf("Testing scalar_ui\n");
            uc = n_randint(state, 0);
            if (l == 0)
                TEST_MUL_SCALAR(status, k, csr, scalar_ui, csr_mat, csr_mat2, dmat, dmat2, uc, ctx)
            else if (l == 1)
                TEST_MUL_SCALAR(status, k, lil, scalar_ui, lil_mat, lil_mat2, dmat, dmat2, uc, ctx)
            else
                TEST_MUL_SCALAR(status, k, coo, scalar_ui, coo_mat, coo_mat2, dmat, dmat2, uc, ctx)
            break;
        case 3:
            //flint_printf("Testing scalar_fmpz\n");
            fmpz_randtest_not_zero(fmpz_c, state, 32);
            if (l == 0)
                TEST_MUL_SCALAR(status, k, csr, scalar_fmpz, csr_mat, csr_mat2, dmat, dmat2, fmpz_c, ctx)
            else if (l == 1)
                TEST_MUL_SCALAR(status, k, lil, scalar_fmpz, lil_mat, lil_mat2, dmat, dmat2, fmpz_c, ctx)
            else
                TEST_MUL_SCALAR(status, k, coo, scalar_fmpz, coo_mat, coo_mat2, dmat, dmat2, fmpz_c, ctx)
            break;
        case 4:
            //flint_printf("Testing scalar_fmpq\n");
            fmpq_randtest_not_zero(fmpq_c, state, 32);
            if (l == 0)
                TEST_MUL_SCALAR(status, k, csr, scalar_fmpq, csr_mat, csr_mat2, dmat, dmat2, fmpq_c, ctx)
            else if (l == 1)
                TEST_MUL_SCALAR(status, k, lil, scalar_fmpq, lil_mat, lil_mat2, dmat, dmat2, fmpq_c, ctx)
            else
                TEST_MUL_SCALAR(status, k, coo, scalar_fmpq, coo_mat, coo_mat2, dmat, dmat2, fmpq_c, ctx)
            break;
        case 5:
            //flint_printf("Testing scalar_2exp_si\n");
            // Scaling by 2^c always done with multiply, even if it is a divide
            c = n_randint(state, 32) + 1;
            if (k == 0 || k == 2)
            {
                if (l == 0)
                { TEST_OP_SCALAR(status, csr, mul, scalar_2exp_si, csr_mat, csr_mat2, dmat, dmat2, c, ctx) }
                else if (l == 1)
                { TEST_OP_SCALAR(status, lil, mul, scalar_2exp_si, lil_mat, lil_mat2, dmat, dmat2, c, ctx) }
                else
                { TEST_OP_SCALAR(status, coo, mul, scalar_2exp_si, coo_mat, coo_mat2, dmat, dmat2, c, ctx) }
            }
            if (k == 1 || k == 2)
            {
                if (l == 0)
                { TEST_OP_SCALAR(status, csr, mul, scalar_2exp_si, csr_mat, csr_mat2, dmat, dmat2, -c, ctx) }
                else if (l == 1)
                { TEST_OP_SCALAR(status, lil, mul, scalar_2exp_si, lil_mat, lil_mat2, dmat, dmat2, -c, ctx) }
                else
                { TEST_OP_SCALAR(status, coo, mul, scalar_2exp_si, coo_mat, coo_mat2, dmat, dmat2, -c, ctx) }
            }
            break;
        }
        
        // If any operation not allowed, just skip test
        if (status == GR_UNABLE || status == GR_DOMAIN) // TODO: FIXME
        {
            status = GR_SUCCESS;
            continue;
        }
        //gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
        //gr_sparse_vec_print_nz(vec2, ctx); flint_printf("\n");
 
        if (l == 0)
            status |= gr_mat_set_csr_mat(dmat, csr_mat2, ctx);
        else if (l == 1)
            status |= gr_mat_set_lil_mat(dmat, lil_mat2, ctx);
        else
            status |= gr_mat_set_coo_mat(dmat, coo_mat2, ctx);
        eq = gr_mat_equal(dmat, dmat2, ctx);
        if (eq == T_FALSE || status != GR_SUCCESS)
        {
            flint_printf(
                "j = %d, k = %d, equal = %d, status = %d\n",
                j, k, eq, status
            );
            gr_ctx_println(ctx);
            flint_printf("\n\ndmat: "); status |= gr_mat_print(dmat, ctx); flint_printf("\n");
            flint_printf("\n\ndmat2: "); status |= gr_mat_print(dmat2, ctx); flint_printf("\n");
            if (l == 0)
            {
                status |= gr_csr_mat_set_mat(csr_mat, dmat2, ctx);
                status |= gr_csr_mat_print_nz(csr_mat, ctx); flint_printf("\n");
                status |= gr_csr_mat_print_nz(csr_mat2, ctx); flint_printf("\n");
            }
            else if (l == 1)
            {
                status |= gr_lil_mat_set_mat(lil_mat, dmat2, ctx);
                status |= gr_lil_mat_print_nz(lil_mat, ctx); flint_printf("\n");
                status |= gr_lil_mat_print_nz(lil_mat2, ctx); flint_printf("\n");
            }
            else
            {
                status |= gr_coo_mat_set_mat(coo_mat, dmat2, ctx);
                status |= gr_coo_mat_print_nz(coo_mat, ctx); flint_printf("\n");
                status |= gr_coo_mat_print_nz(coo_mat2, ctx); flint_printf("\n");
            }
            return GR_TEST_FAIL;
        }
    }

    fmpz_clear(fmpz_c);
    fmpq_clear(fmpq_c);
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

TEST_FUNCTION_START(gr_sparse_mat_arith, state)
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
        CHECK_TEST(test_neg(state, ctx), "Sparse matrix negation");
        CHECK_TEST(test_add_sub(state, ctx), "Sparse matrix addition, subtraction, and multiplication");
        CHECK_TEST(test_accum_mul_scalar(state, ctx), "Sparse matrix scalar addmul and submul");
        CHECK_TEST(test_mul_div_scalar(state, ctx), "Sparse matrix scalar multiplication and division");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
