/*
    Copyright (C) 2023 Vincent Neiger
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"
#include "gr_mat.h"

TEST_GR_FUNCTION_START(gr_mat_permute_rows, state, count_success, count_domain, count_unable)
{
    slong iter;
    
    
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong m, n, i, j;
        gr_ctx_t ctx;
        gr_mat_t mat, mat2, mat3;
        slong *perm_act, *perm_act_inv, *perm, *perm_store, *perm_store_copy;

        slong status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        m = n_randint(state, 6);
        n = n_randint(state, 6);

        gr_mat_init(mat, m, n, ctx);
        status |= gr_mat_randtest(mat, state, ctx);

        perm_act = _perm_init(m);
        _perm_randtest(perm_act, m, state);
        perm_store = _perm_init(m);
        _perm_randtest(perm_store, m, state);
        perm = _perm_init(m); /* copy of perm_store for testing purpose */
        _perm_set(perm, perm_store, m);

        status |= gr_mat_init_set(mat2, mat, ctx);
        status |= gr_mat_permute_rows(mat2, perm_store, perm_act, ctx);

        for (i = 0; i < m; i++)
        {
            if (perm_store && perm_store[i] != perm[perm_act[i]])
            {
                flint_printf("FAIL:\n");
                flint_printf("auxiliary permutation not correctly permuted by perm_act\n");
                flint_printf("m = %wd\n", m);
                flint_printf("input permutation: %{slong*}\n", perm, m);
                flint_printf("acting permutation: %{slong*}\n", perm_act, m);
                flint_printf("resulting permutation: %{slong*}\n", perm_store, m);
                flint_abort();
            }

            for (j = 0; j < n; j++)
            {
                if (gr_equal(GR_MAT_ENTRY(mat2, i, j, ctx->sizeof_elem), GR_MAT_ENTRY(mat, perm_act[i], j, ctx->sizeof_elem), ctx) == T_FALSE)
                {
                    flint_printf("FAIL (2):\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    flint_printf("matrix not correctly row-permuted by perm_act\n");
                    flint_printf("first matrix = "); gr_mat_print(mat, ctx); flint_printf("\n");
                    flint_printf("acting permutation: %{slong*}\n", perm_act, m);
                    flint_printf("second matrix = "); gr_mat_print(mat2, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        // Now permute back using the inverse of perm_act (in two different ways) to verify that the operation is reversible
        status |= gr_mat_init_set(mat3, mat2, ctx);
        perm_act_inv = _perm_init(m);
        _perm_inv(perm_act_inv, perm_act, m);
        perm_store_copy = _perm_init(m);
        _perm_set(perm_store_copy, perm_store, m);
        
        status |= gr_mat_permute_rows_inv(mat2, perm_store, perm_act, ctx);

        if (!_perm_equal(perm_store, perm, m))
        {
            flint_printf("FAIL (3):\n");
            flint_printf("auxiliary permutation not correctly permuted back by gr_mat_permute_rows_inv\n");
            flint_abort();
        }

        if (gr_mat_equal(mat, mat2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (4):\n");
            flint_printf("matrix not correctly row-permuted back by gr_mat_permute_rows_inv\n");
            flint_abort();
        }

        _perm_set(perm_store, perm_store_copy, m);
        status |= gr_mat_permute_rows(mat3, perm_store, perm_act_inv, ctx);

        if (!_perm_equal(perm_store, perm, m))
        {
            flint_printf("FAIL (5):\n");
            flint_printf("auxiliary permutation not correctly permuted back by inverse of perm_act\n");
            flint_abort();
        }

        if (gr_mat_equal(mat, mat2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (6):\n");
            flint_printf("matrix not correctly row-permuted back by inverse of perm_act\n");
            flint_abort();
        }

        // Act once again with the original permutation, but this time with perm_store set to one
        _perm_one(perm_store, m);
        status |= gr_mat_permute_rows(mat2, perm_store, perm_act, ctx);
        status |= gr_mat_permute_rows_inv(mat2, perm_store, perm_act, ctx);
        if (!_perm_is_one(perm_store, m))
        {
            flint_printf("FAIL (7):\n");
            flint_printf("auxiliary permutation not correctly permuted back by its inverse\n");
            flint_abort();
        }

        if (gr_mat_equal(mat, mat2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (8):\n");
            flint_printf("matrix not correctly row-permuted back by inverse of perm_store\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(mat, ctx);
        gr_mat_clear(mat2, ctx);
        gr_mat_clear(mat3, ctx);
        gr_ctx_clear(ctx);
        _perm_clear(perm_act);
        _perm_clear(perm_act_inv);
        _perm_clear(perm);
        _perm_clear(perm_store);
        _perm_clear(perm_store_copy);
    }

    // consistency with gr_mat_swap_rows
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong m = 3;
        slong n = 3;
        gr_ctx_t ctx;
        gr_mat_t mat, mat2;
        slong * perm;

        slong status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        gr_mat_init(mat, m, n, ctx);
        status |= gr_mat_randtest(mat, state, ctx);

        status |= gr_mat_init_set(mat2, mat, ctx);

        perm = _perm_init(m);
        _perm_one(perm, m);

        status |= gr_mat_swap_rows(mat2, perm, 0, 1, ctx);
        status |= gr_mat_swap_rows(mat2, perm, 1, 2, ctx);

        status |= gr_mat_permute_rows_inv(mat2, perm, perm, ctx);

        if (!_perm_is_one(perm, m))
        {
            flint_printf("FAIL (B 1):\n");
            flint_printf("auxiliary permutation not correctly permuted back by its inverse\n");
            flint_abort();
        }

        if (gr_mat_equal(mat, mat2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (B 2):\n");
            flint_printf("matrix not correctly row-permuted back by inverse of perm\n");
            flint_abort();
        }
        
        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(mat, ctx);
        gr_mat_clear(mat2, ctx);
        gr_ctx_clear(ctx);
        _perm_clear(perm);
    }


    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
