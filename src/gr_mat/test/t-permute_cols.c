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

TEST_GR_FUNCTION_START(gr_mat_permute_cols, state, count_success, count_domain, count_unable)
{
    slong iter;
    
    
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong m, n, i, j;
        gr_ctx_t ctx;
        gr_mat_t mat, mat2;
        slong * perm_act;
        slong * perm_store;
        slong * perm;

        slong status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        m = n_randint(state, 6);
        n = n_randint(state, 6);

        gr_mat_init(mat, m, n, ctx);
        status |= gr_mat_randtest(mat, state, ctx);

        perm_act = _perm_init(n);
        _perm_randtest(perm_act, n, state);
        perm_store = _perm_init(n);
        _perm_randtest(perm_store, n, state);
        perm = _perm_init(n); /* copy of perm_store for testing purpose */
        _perm_set(perm, perm_store, n);

        status |= gr_mat_init_set(mat2, mat, ctx);
        status |= gr_mat_permute_cols(mat2, perm_store, perm_act, ctx);

        for (j = 0; j < n; j++)
        {
            if (perm_store && perm_store[j] != perm[perm_act[j]])
            {
                flint_printf("FAIL:\n");
                flint_printf("auxiliary permutation not correctly permuted by perm_act\n");
                flint_printf("n = %wd\n", n);
                flint_printf("input permutation: %{slong*}\n", perm, n);
                flint_printf("acting permutation: %{slong*}\n", perm_act, n);
                flint_printf("resulting permutation: %{slong*}\n", perm_store, n);
                flint_abort();
            }

            for (i = 0; i < m; i++)
            {
                if (gr_equal(GR_MAT_ENTRY(mat2, i, j, ctx->sizeof_elem), GR_MAT_ENTRY(mat, i, perm_act[j], ctx->sizeof_elem), ctx) == T_FALSE)
                {
                    flint_printf("FAIL (2):\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    flint_printf("matrix not correctly col-permuted by perm_act\n");
                    flint_printf("first matrix = "); gr_mat_print(mat, ctx); flint_printf("\n");
                    flint_printf("acting permutation: %{slong*}\n", perm_act, m);
                    flint_printf("second matrix = "); gr_mat_print(mat2, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        // Now permute back using the inverse of perm_act to verify that the operation is reversible
        _perm_inv(perm_act, perm_act, n);
        status |= gr_mat_permute_cols(mat2, perm_store, perm_act, ctx);
        
        if (!_perm_equal(perm_store, perm, n))
        {
            flint_printf("FAIL (3):\n");
            flint_printf("auxiliary permutation not correctly permuted back by inverse of perm_act\n");
            flint_abort();
        }

        if (gr_mat_equal(mat, mat2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (4):\n");
            flint_printf("matrix not correctly col-permuted back by inverse of perm_act\n");
            flint_abort();
        }

        // Act once again with the original permutation, but this time with perm_store set to one
        _perm_inv(perm_act, perm_act, n);
        _perm_one(perm_store, n);
        status |= gr_mat_permute_cols(mat2, perm_store, perm_act, ctx);
        _perm_inv(perm_act, perm_store, n);
        status |= gr_mat_permute_cols(mat2, perm_store, perm_act, ctx);
        if (!_perm_is_one(perm_store, n))
        {
            flint_printf("FAIL (5):\n");
            flint_printf("auxiliary permutation not correctly permuted back by its inverse\n");
            flint_abort();
        }

        if (gr_mat_equal(mat, mat2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (B 2):\n");
            flint_printf("matrix not correctly col-permuted back by inverse of perm_store\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(mat, ctx);
        gr_mat_clear(mat2, ctx);
        gr_ctx_clear(ctx);
        _perm_clear(perm_act);
        _perm_clear(perm);
        _perm_clear(perm_store);
    }

    // consistency with gr_mat_swap_cols
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong m = 3;
        slong n = 3;
        gr_ctx_t ctx;
        gr_mat_t mat, mat2;
        slong * perm_act;
        slong * perm_store;

        slong status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        gr_mat_init(mat, m, n, ctx);
        status |= gr_mat_randtest(mat, state, ctx);

        status |= gr_mat_init_set(mat2, mat, ctx);

        perm_store = _perm_init(n);
        _perm_one(perm_store, n);

        status |= gr_mat_swap_cols(mat2, perm_store, 0, 1, ctx);
        status |= gr_mat_swap_cols(mat2, perm_store, 1, 2, ctx);

        perm_act = _perm_init(n);
        _perm_inv(perm_act, perm_store, n);

        status |= gr_mat_permute_cols(mat2, perm_store, perm_act, ctx);

        if (!_perm_is_one(perm_store, n))
        {
            flint_printf("FAIL (B 1):\n");
            flint_printf("auxiliary permutation not correctly permuted back by its inverse\n");
            flint_abort();
        }

        if (gr_mat_equal(mat, mat2, ctx) == T_FALSE)
        {
            flint_printf("FAIL (B 2):\n");
            flint_printf("matrix not correctly col-permuted back by inverse of perm_store\n");
            flint_abort();
        }
        
        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(mat, ctx);
        gr_mat_clear(mat2, ctx);
        gr_ctx_clear(ctx);
        _perm_clear(perm_act);
        _perm_clear(perm_store);
    }


    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
