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
        gr_mat_t mat, matt;
        slong * perm_act;
        slong * perm_store;
        slong * perm;

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

        status |= gr_mat_init_set(matt, mat, ctx);
        status |= gr_mat_permute_rows(matt, perm_store, perm_act, ctx);

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
                if (gr_equal(GR_MAT_ENTRY(matt, i, j, ctx->sizeof_elem), GR_MAT_ENTRY(mat, perm_act[i], j, ctx->sizeof_elem), ctx) == T_FALSE)
                {
                    flint_printf("FAIL (2):\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    flint_printf("matrix not correctly row-permuted by perm_act\n");
                    flint_printf("first matrix = "); gr_mat_print(mat, ctx); flint_printf("\n");
                    flint_printf("acting permutation: %{slong*}\n", perm_act, m);
                    flint_printf("second matrix = "); gr_mat_print(matt, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(mat, ctx);
        gr_mat_clear(matt, ctx);
        gr_ctx_clear(ctx);
        _perm_clear(perm_act);
        _perm_clear(perm);
        _perm_clear(perm_store);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
