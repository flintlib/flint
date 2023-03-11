/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, mod, rep;
    FLINT_TEST_INIT(state);


    flint_printf("permute_rows....");
    fflush(stdout);

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
                flint_printf("FAIL: auxiliary permutation not correctly permuted by perm_act\n");
                flint_printf("input permutation:\n");
                _perm_print(perm, m);
                flint_printf("acting permutation:\n");
                _perm_print(perm_act, m);
                flint_printf("resulting permutation:\n");
                _perm_print(perm_store, m);
                fflush(stdout);
                flint_abort();
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

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
