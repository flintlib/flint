/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_scalar_mul, state)
{
    slong m, n, mod, rep;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C, D;
        ulong c;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_not_zero(state);

        c = n_randint(state, mod);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(C, m, n, mod);
        nmod_mat_init(D, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_scalar_mul(C, A, c);
        nmod_mat_scalar_mul(D, A, nmod_sub(c, nmod_set_ui(1, A->mod), A->mod));

        /* c*A - (c-1)*A == A */
        nmod_mat_sub(D, C, D);

        if (!nmod_mat_equal(A, D))
            TEST_FUNCTION_FAIL("");

        /* Aliasing */
        nmod_mat_scalar_mul(C, A, c);
        nmod_mat_scalar_mul(A, A, c);

        if (!nmod_mat_equal(A, C))
            TEST_FUNCTION_FAIL("");

        /* entry by entry, into a window of a larger matrix (rows that do
           not follow each other) and from one */
        {
            nmod_mat_t big, W;
            slong i, j, r0 = n_randint(state, 3), c0 = n_randint(state, 4);

            nmod_mat_init(big, m + 2, n + 3, mod);
            nmod_mat_randtest(big, state);
            nmod_mat_window_init(W, big, r0, c0, r0 + m, c0 + n);
            nmod_mat_scalar_mul(W, B, c);
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    if (nmod_mat_entry(W, i, j) != nmod_mul(c, nmod_mat_entry(B, i, j), B->mod))
                        TEST_FUNCTION_FAIL("window: i = %wd, j = %wd\n", i, j);
            nmod_mat_set(D, W);
            nmod_mat_scalar_mul(C, W, c);
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    if (nmod_mat_entry(C, i, j) != nmod_mul(c, nmod_mat_entry(D, i, j), B->mod))
                        TEST_FUNCTION_FAIL("from a window: i = %wd, j = %wd\n", i, j);
            nmod_mat_window_clear(W);
            nmod_mat_clear(big);
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
