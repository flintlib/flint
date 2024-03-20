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

int test_init_csr(gr_ctx_t ctx)
{
    gr_csr_mat_t mat;
    gr_csr_mat_init(mat, 5, 6, ctx);
    if (mat->r != 5 || mat->c != 6 || mat->alloc != 0 || mat->nnz != 0 || mat->rows == NULL || mat->cols != NULL || mat->nzs != NULL)
    {
        flint_printf("Failed init CSR test\n");
        return GR_TEST_FAIL;
    }
    gr_csr_mat_clear(mat, ctx);
    return GR_SUCCESS;
}

int test_init_lil(gr_ctx_t ctx)
{
    gr_lil_mat_t mat;
    gr_lil_mat_init(mat, 5, 6, ctx);
    if (mat->r != 5 || mat->c != 6 || mat->nnz != 0 || mat->rows == NULL)
    {
        flint_printf("Failed init LIL test\n");
        return GR_TEST_FAIL;
    }
    gr_lil_mat_clear(mat, ctx);
    return GR_SUCCESS;
}

int test_init_coo(gr_ctx_t ctx)
{
    gr_coo_mat_t mat;
    gr_coo_mat_init(mat, 5, 6, ctx);
    if (mat->r != 5 || mat->c != 6 || mat->alloc != 0 || mat->nnz != 0 || mat->rows != NULL || mat->cols != NULL || mat->nzs != NULL || mat->is_canonical != T_TRUE)
    {
        flint_printf("Failed init COO test\n");
        return GR_TEST_FAIL;
    }
    gr_coo_mat_clear(mat, ctx);
    return GR_SUCCESS;
}

int test_init_from_entries_canonical(flint_rand_t state, gr_ctx_t ctx)
{
    // TODO: randomize length and nonzero cols
    slong i;
    gr_coo_mat_t mat;
    gr_ptr entries;
    gr_ptr temp;
    int status = GR_SUCCESS;
    truth_t eq;
    slong sz = ctx->sizeof_elem;
    slong r = 5;
    slong c = 10;
    slong N = 15;
    ulong rows[25] = {0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4};
    ulong cols[25] =  {0, 2, 3, 1, 4, 6, 8, 9, 5, 7, 0, 1, 3, 5, 6, 8, 2, 4, 7, 9};

    //flint_printf("Running init test\n");
    GR_TMP_INIT(temp, ctx);
    GR_TMP_INIT_VEC(entries, N, ctx);
    for (i = 0; i < N; ++i)
        status |= gr_randtest_not_zero(GR_ENTRY(entries, i, sz), state, ctx);
    if (status != GR_SUCCESS)
    {
        flint_printf("Failed to make random numbers!\n");
        return GR_SUCCESS; // Not my fault!
    }

    gr_coo_mat_init(mat, r, c, ctx);
    status |= gr_coo_mat_from_entries(mat, rows, cols, entries, N, T_TRUE, ctx);
    
    // Check parameters
    if (status != GR_SUCCESS || mat->r != r || mat->c != c || mat->alloc != 15 || mat->nnz != 15)
    {
        flint_printf("Bad params! %ld %ld %ld %ld\n", mat->r, mat->c, mat->alloc, mat->nnz);
        return GR_TEST_FAIL;
    }

    // Check indices and entries
    for (i = 0; i < N; ++i)
    {
        if (*gr_coo_mat_row_ptr(mat, i) != rows[i])
        {
            flint_printf("Bad row index!\n");
            return GR_TEST_FAIL;
        }
        if (*gr_coo_mat_col_ptr(mat, i) != cols[i])
        {
            flint_printf("Bad column index!\n");
            return GR_TEST_FAIL;
        }
        eq = gr_equal(gr_coo_mat_entry_ptr(mat, i, ctx), GR_ENTRY(entries, i, sz), ctx);
        if (eq == T_FALSE)
        {
            flint_printf("Bad elements!\n");
            return GR_TEST_FAIL;
        }
    }
    gr_coo_mat_clear(mat, ctx);
    GR_TMP_CLEAR_VEC(entries, N, ctx);
    return status;
}

int test_init_from_entries_internal(ulong *rows, ulong *cols, gr_srcptr entries, slong r, slong c, slong num, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, j, k;
    gr_coo_mat_t mat;
    gr_ptr temp, temp2, temp3;

    ///flint_printf("\n\n\nRunning internal test\n");
    GR_TMP_INIT2(temp, temp2, ctx);
    gr_coo_mat_init(mat, r, c, ctx);
    status |= gr_coo_mat_from_entries(mat, rows, cols, entries, num, T_FALSE, ctx);
    if (status != GR_SUCCESS)
        return GR_TEST_FAIL;

    gr_coo_mat_canonicalize(mat, ctx);
    if (mat->is_canonical == T_FALSE)
        return GR_TEST_FAIL;

    // Check every entry (including the zeroes)
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            // Compute the expected value of the entry
            status |= gr_zero(temp, ctx);
            for (k = 0; k < num; k++)
                if (rows[k] == i && cols[k] == j)
                    status |= gr_add(temp, temp, GR_ENTRY(entries, k, sz), ctx);

            status |= gr_coo_mat_get_entry(temp2, mat, i, j, ctx);
            temp3 = gr_coo_mat_find_entry(mat, i, j, ctx);
            if (
                gr_equal(temp, temp2, ctx) == T_FALSE || 
                (temp3 == NULL && gr_is_zero(temp, ctx) == T_FALSE) ||
                (temp3 != NULL && gr_is_zero(temp3, ctx) == T_TRUE) ||
                (temp3 != NULL && gr_equal(temp, temp3, ctx) == T_FALSE)
            )
                {
                    flint_printf("Failed on %d!\n", i);
                    gr_println(temp, ctx);
                    gr_println(temp2, ctx);
                    flint_printf("%p\n", temp3);
                    if (temp3 != NULL)
                        gr_println(temp3, ctx);
                    return GR_TEST_FAIL;
                }
        }
    }
    GR_TMP_CLEAR2(temp, temp2, ctx);
    gr_coo_mat_clear(mat, ctx);
    return status;
}

int test_init_from_entries(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, j;
    int status = GR_SUCCESS;
    slong r = 5;
    slong c = 10;
    slong N = 20;
    ulong rows[20] = {0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4};
    ulong cols[20] =  {0, 2, 3, 1, 4, 6, 8, 9, 5, 7, 0, 1, 3, 5, 6, 8, 2, 4, 7, 9};
    gr_ptr entries;
    GR_TMP_INIT_VEC(entries, N, ctx);
    status |= _gr_vec_randtest(entries, state, N, ctx);

    // Randomly permute the rows and columns
    for (i = 0; i < N; ++i)
    {
        j = n_randint(state, N);
        FLINT_SWAP(ulong, rows[i], rows[j]);
        j = n_randint(state, N);
        FLINT_SWAP(ulong, cols[i], cols[j]);
    }
    status |= test_init_from_entries_internal(rows, cols, entries, r, c, N, ctx);

    GR_TMP_CLEAR_VEC(entries, N, ctx);
    return status;
}

int test_init_from_entries_adversarial(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong r = 5;
    slong c = 10;
    slong N = 20;
    ulong rows[20] = {1,  0, 3, 4, 2, 1, 4, 1, 3, 2,  3,  3, 4, 0, 0,  1, 2, 4,  3, 0};
    ulong cols[20] = {0,  2, 3, 1, 4, 8, 8, 0, 5, 7,  0,  3, 3, 5, 6,  8, 4, 2,  5, 9};
    slong ents[20] = {-1, 3, 1, 2, 0, 6, 4, 1, 3, 1, -4, -1, 2, 9, 1, -6, 4, 0, -2, 4};
    
    gr_ptr entries;
    GR_TMP_INIT_VEC(entries, N, ctx);
    
    for (i = 0; i < N; i++)
        status |= gr_set_si(GR_ENTRY(entries, i, sz), ents[i], ctx);
    status |= test_init_from_entries_internal(rows, cols, entries, r, c, N, ctx);

    GR_TMP_CLEAR_VEC(entries, N, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_mat_init, state)
{   
    int i;
    gr_ctx_t ctx;
    for (i = 0; i < 16; ++i)
    {
        while(1)
        {
            //gr_ctx_init_fmpz(ctx);
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        //gr_ctx_println(ctx);
        CHECK_TEST(test_init_csr(ctx), "Init CSR matrix")
        CHECK_TEST(test_init_lil(ctx), "Init LIL matrix")
        CHECK_TEST(test_init_coo(ctx), "Init COO matrix")
        CHECK_TEST(test_init_from_entries_canonical(state, ctx), "Init from entries in canonical form");
        CHECK_TEST(test_init_from_entries(state, ctx), "Init from entries in noncanonical order");
        CHECK_TEST(test_init_from_entries_adversarial(state, ctx), "Init from entries in adversarial order");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
