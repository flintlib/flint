/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_permute_rows, state)
{
    slong m, n, mod, rep;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t mat, matt;
        slong * perm_act;
				slong * perm_store;
				slong * perm;
        slong i, j;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(mat, m, n, mod);
        nmod_mat_randtest(mat, state);

				perm_act = _perm_init(m);
        _perm_randtest(perm_act, m, state);
				perm_store = _perm_init(m);
        _perm_randtest(perm_store, m, state);
        /* perm = copy of perm_store for testing purpose */
				perm = _perm_init(m);
        for (i = 0; i < m; i++)
            perm[i] = perm_store[i];

        nmod_mat_init_set(matt, mat);
        nmod_mat_permute_rows(matt, perm_act, perm_store);

        for (i = 0; i < m; i++)
        {
            if (perm_store && perm_store[i] != perm[perm_act[i]])
            {
                flint_throw(FLINT_TEST_FAIL,
                        "auxiliary permutation not correctly permuted by perm_act\n"
                        "m = %wd\n"
                        "input permutation: %{slong*}\n"
                        "acting permutation: %{slong*}\n"
                        "resulting permutation: %{slong*}\n",
                        m,
                        perm, m,
                        perm_act, m,
                        perm_store, m);
            }
            for (j = 0; j < n; j++)
            {
                if (nmod_mat_entry(matt, i, j) != nmod_mat_entry(mat, perm_act[i], j))
                {
                    flint_printf("FAIL: matrix not correctly row-permuted by perm_act\n");
                    flint_printf("first matrix:\n");
                    nmod_mat_print_pretty(mat);
                    flint_printf("second matrix:\n");
                    nmod_mat_print_pretty(matt);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        nmod_mat_clear(mat);
        nmod_mat_clear(matt);
        _perm_clear(perm_act);
        _perm_clear(perm);
        _perm_clear(perm_store);
    }

    TEST_FUNCTION_END(state);
}
