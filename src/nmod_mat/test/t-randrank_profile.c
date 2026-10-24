/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "rank_profile_helpers.h"

TEST_FUNCTION_START(nmod_mat_randrank_profile, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        nmod_mat_t A, B, W;
        slong m, n, rank, r2, k, i, j;
        slong * mrp_r, * mrp_c, * mrp_r2, * mrp_c2, * sh_r, * sh_c;
        ulong p;
        int window;

        m = n_randint(state, 12);
        n = n_randint(state, 12);
        if (n_randint(state, 50) == 0)
        {
            m = n_randint(state, 30);
            n = n_randint(state, 30);
        }
        p = _rp_rand_prime(state);
        rank = n_randint(state, FLINT_MIN(m, n) + 1);
        window = n_randint(state, 4) == 0;

        mrp_r = flint_malloc((rank + 1) * sizeof(slong));
        mrp_c = flint_malloc((rank + 1) * sizeof(slong));
        mrp_r2 = flint_malloc((FLINT_MIN(m, n) + 1) * sizeof(slong));
        mrp_c2 = flint_malloc((FLINT_MIN(m, n) + 1) * sizeof(slong));
        sh_r = flint_malloc((rank + 1) * sizeof(slong));
        sh_c = flint_malloc((rank + 1) * sizeof(slong));
        _rp_rand_positions(mrp_r, mrp_c, m, n, rank, state);

        /* the positions given in any order (mrp_r, mrp_c are sorted by
           column, as the rank profile matrix computed below) */
        for (k = 0; k < rank; k++)
        {
            sh_r[k] = mrp_r[k];
            sh_c[k] = mrp_c[k];
        }
        for (k = rank - 1; k > 0; k--)
        {
            i = n_randint(state, k + 1);
            FLINT_SWAP(slong, sh_r[k], sh_r[i]);
            FLINT_SWAP(slong, sh_c[k], sh_c[i]);
        }

        nmod_mat_init(A, m, n, p);
        nmod_mat_randtest(A, state);

        if (window)
        {
            /* generate in a window of a larger matrix B, whose other
               entries must be preserved */
            nmod_mat_t B0;
            nmod_mat_init(B, m + 2, n + 3, p);
            nmod_mat_randtest(B, state);
            nmod_mat_init_set(B0, B);
            nmod_mat_window_init(W, B, 1, 2, m + 1, n + 2);
            nmod_mat_randrank_profile(W, state, sh_r, sh_c, rank);
            nmod_mat_set(A, W);
            nmod_mat_window_clear(W);
            for (i = 0; i < m + 2; i++)
                for (j = 0; j < n + 3; j++)
                    if ((i < 1 || i >= m + 1 || j < 2 || j >= n + 2)
                        && nmod_mat_entry(B, i, j) != nmod_mat_entry(B0, i, j))
                        TEST_FUNCTION_FAIL("window: entry (%wd, %wd) modified\n", i, j);
            nmod_mat_clear(B);
            nmod_mat_clear(B0);
        }
        else
        {
            nmod_mat_randrank_profile(A, state, sh_r, sh_c, rank);
        }

        r2 = _rp_rank_profile_naive(mrp_r2, mrp_c2, A);

        if (r2 != rank)
            TEST_FUNCTION_FAIL(
                    "rank: %wd, expected %wd\n"
                    "A = %{nmod_mat}\n", r2, rank, A);

        for (k = 0; k < rank; k++)
            if (mrp_r[k] != mrp_r2[k] || mrp_c[k] != mrp_c2[k])
                TEST_FUNCTION_FAIL(
                        "rank profile matrix, pivot %wd: (%wd, %wd), expected (%wd, %wd)\n"
                        "A = %{nmod_mat}\n",
                        k, mrp_r2[k], mrp_c2[k], mrp_r[k], mrp_c[k], A);

        nmod_mat_clear(A);
        flint_free(mrp_r);
        flint_free(mrp_c);
        flint_free(mrp_r2);
        flint_free(mrp_c2);
        flint_free(sh_r);
        flint_free(sh_c);
    }

    TEST_FUNCTION_END(state);
}
