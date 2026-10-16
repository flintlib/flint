/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"

TEST_FUNCTION_START(nmod_mat_mul_strassen, state)
{
    slong i;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        ulong mod = n_randtest_not_zero(state);

        slong m, k, n;

        m = n_randint(state, 400);
        k = n_randint(state, 400);
        n = n_randint(state, 400);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, k, mod);
        nmod_mat_init(C, m, k, mod);
        nmod_mat_init(D, m, k, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_mul_classical(C, A, B);
        nmod_mat_mul_strassen(D, A, B);

        if (!nmod_mat_equal(C, D))
            TEST_FUNCTION_FAIL(
                    "Results not equal\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "C = %{nmod_mat}\n"
                    "D = %{nmod_mat}\n",
                    A, B, C, D);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    /* several levels with small cutoffs: all parities of the dimensions
       at every level (virtual padding), operands and result as windows */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D, PA, PB, PC;
        slong m, k, n, cutoff, r, c;
        ulong mod;

        m = 5 + n_randint(state, 150);
        k = 5 + n_randint(state, 150);
        n = 5 + n_randint(state, 150);
        cutoff = 5 + n_randint(state, 30);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(PA, m + 2, k + 3, mod);
        nmod_mat_init(PB, k + 1, n + 2, mod);
        nmod_mat_init(PC, m + 3, n + 1, mod);
        nmod_mat_randtest(PA, state);
        nmod_mat_randtest(PB, state);
        nmod_mat_randtest(PC, state);
        nmod_mat_window_init(A, PA, 1, 2, m + 1, k + 2);
        nmod_mat_window_init(B, PB, 1, 1, k + 1, n + 1);
        nmod_mat_window_init(C, PC, 2, 0, m + 2, n);
        nmod_mat_init(D, m + 3, n + 1, mod);
        nmod_mat_set(D, PC);

        _nmod_mat_mul_strassen_cutoff(C, A, B, cutoff);

        /* D: the expected parent of C */
        {
            nmod_mat_t W;
            nmod_mat_window_init(W, D, 2, 0, m + 2, n);
            nmod_mat_mul_classical(W, A, B);
            nmod_mat_window_clear(W);
        }

        for (r = 0; r < PC->r; r++)
            for (c = 0; c < PC->c; c++)
                if (nmod_mat_entry(PC, r, c) != nmod_mat_entry(D, r, c))
                    TEST_FUNCTION_FAIL("recursive, m: %wd, k: %wd, n: %wd, "
                                       "mod: %wu, cutoff: %wd, entry %wd, %wd\n",
                                       m, k, n, mod, cutoff, r, c);

        nmod_mat_window_clear(A);
        nmod_mat_window_clear(B);
        nmod_mat_window_clear(C);
        nmod_mat_clear(D);
        nmod_mat_clear(PA);
        nmod_mat_clear(PB);
        nmod_mat_clear(PC);
    }

    TEST_FUNCTION_END(state);
}
