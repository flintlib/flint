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

int test_init(gr_ctx_t ctx)
{
    gr_sparse_vec_t vec;
    gr_sparse_vec_init(vec, 5, ctx);
    if (vec->length != 5 || vec->alloc != 0 || vec->nnz != 0 || vec->inds != NULL || vec->nzs != NULL)
        return GR_TEST_FAIL;
    gr_sparse_vec_clear(vec, ctx);
    return GR_SUCCESS;
}

int
test_init_from_entries_canonical(flint_rand_t state, gr_ctx_t ctx)
{
    // TODO: randomize length and nonzero cols
    slong i;
    gr_sparse_vec_t vec;
    gr_ptr entries;
    gr_ptr temp;
    int status = GR_SUCCESS;
    truth_t eq;
    slong sz = ctx->sizeof_elem;
    slong N = 5;
    slong len = 10;
    ulong inds[5] = {0, 2, 3, 6, 9};

    //flint_printf("Running init test\n");
    GR_TMP_INIT(temp, ctx);
    GR_TMP_INIT_VEC(entries, N, ctx);
    for (i = 0; i < N; ++i)
        status |= gr_randtest_not_zero(GR_ENTRY(entries, i, sz), state, ctx);
    if (status != GR_SUCCESS)
    {
        flint_printf("Failed to make random numbers!\n");
        return GR_TEST_FAIL; // Not my fault!
    }

    gr_sparse_vec_init(vec, len, ctx);
    status |= gr_sparse_vec_from_entries(vec, inds, entries, N, 1, ctx);
    
    // Check parameters
    if (status != GR_SUCCESS || vec->length != len || vec->alloc != 5 || vec->nnz != 5)
    {
        flint_printf("Bad params! %ld %ld %ld\n", vec->length, vec->alloc, vec->nnz);
        return GR_TEST_FAIL;
    }

    // Check indices and entries
    for (i = 0; i < N; ++i)
    {
        if (*gr_sparse_vec_ind_ptr(vec, i, ctx) != inds[i])
        {
            flint_printf("Bad indices!\n");
            return GR_TEST_FAIL;
        }
        eq = gr_equal(gr_sparse_vec_entry_ptr(vec, i, ctx), GR_ENTRY(entries, i, sz), ctx);
        if (eq == T_FALSE)
        {
            flint_printf("Bad elements!\n");
            return GR_TEST_FAIL;
        }
    }
    gr_sparse_vec_clear(vec, ctx);
    GR_TMP_CLEAR_VEC(entries, N, ctx);
    return status;
}

int
test_init_from_entries_internal(ulong *inds, gr_srcptr entries, slong len, slong num, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i,j;
    gr_sparse_vec_t vec;
    gr_ptr temp, temp2, temp3;

    GR_TMP_INIT2(temp, temp2, ctx);
    gr_sparse_vec_init(vec, len, ctx);
    //flint_printf("entries: "); status |= _gr_vec_print(entries, num, ctx); flint_printf("\n");
    status |= gr_sparse_vec_from_entries(vec, inds, entries, num, T_FALSE, ctx);
    //flint_printf("vec: "); status |= gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
    if (status != GR_SUCCESS)
        return GR_TEST_FAIL;

    // Check every entry (including the zeroes)
    for (i = 0; i < len; i++)
    {
        // Compute the expected value of the entry
        status |= gr_zero(temp, ctx);
        for (j = 0; j < num; j++)
            if (inds[j] == i)
                status |= gr_add(temp, temp, GR_ENTRY(entries, j, sz), ctx);

        status |= gr_sparse_vec_get_entry(temp2, vec, i, ctx);
        temp3 = gr_sparse_vec_find_entry(vec, i, ctx);
        if (
            gr_equal(temp, temp2, ctx) == T_FALSE || 
            (temp3 == NULL && gr_is_zero(temp, ctx) == T_FALSE) ||
            (temp3 != NULL && gr_is_zero(temp3, ctx) == T_TRUE) ||
            (temp3 != NULL && gr_equal(temp, temp3, ctx) == T_FALSE)
        )
            {
                flint_printf("Failed on %d!\n", i);
                gr_ctx_println(ctx);
                gr_println(temp, ctx);
                gr_println(temp2, ctx);
                if (temp3 != NULL)
                    gr_println(temp3, ctx);
                return GR_TEST_FAIL;
            }
    }
    GR_TMP_CLEAR(temp, ctx);
    GR_TMP_CLEAR(temp2, ctx);
    gr_sparse_vec_clear(vec, ctx);
    return status;
}

int
test_init_from_entries(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong N = 5;
    ulong inds[5] = {8, 4, 3, 8, 1};
    gr_ptr entries;

    GR_TMP_INIT_VEC(entries, N, ctx);

    status |= _gr_vec_randtest(entries, state, N, ctx);
    status |= test_init_from_entries_internal(inds, entries, 2*N, N, ctx);

    /* Next test against some adversarial entries */
    slong entries_si[5] = {5, 0, 2, -5, 1};
    for (i = 0; i < N; i++)
        status |= gr_set_si(GR_ENTRY(entries, i, sz), entries_si[i], ctx);
    status |= test_init_from_entries_internal(inds, entries, 2*N, N, ctx);

    GR_TMP_CLEAR_VEC(entries, N, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_vec_init, state)
{   
    int i;
    gr_ctx_t ctx;
    for (i = 0; i < 16; ++i)
    {
        while(1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        //gr_ctx_print(ctx); flint_printf("\n");
        CHECK_TEST(test_init(ctx), "Init");
        CHECK_TEST(test_init_from_entries_canonical(state, ctx), "Init from entries in canonical form");
        CHECK_TEST(test_init_from_entries(state, ctx), "Init from entries in noncanonical form");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
