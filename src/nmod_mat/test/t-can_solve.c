/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_can_solve, state)
{
    nmod_mat_t A, X, X2, B, AX;
    slong i, k, m, n;
    mp_limb_t mod;
    int solved;

    /* test random systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        k = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, k, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(X, k, n, mod);
        nmod_mat_init(AX, m, n, mod);

        nmod_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1));
        nmod_mat_randtest(B, state);

        /* Dense */
        if (n_randint(state, 2))
            nmod_mat_randops(A, state, 1+n_randint(state, 1+m*m));

        solved = nmod_mat_can_solve(X, A, B);

        nmod_mat_mul(AX, A, X);

        if (solved && !nmod_mat_equal(AX, B))
            TEST_FUNCTION_FAIL(
                    "AX != B\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "X = %{nmod_mat}\n"
                    "AX = %{nmod_mat}\n",
                    A, B, X, AX);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(X);
        nmod_mat_clear(AX);
    }

    /* test random solvable systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        k = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, k, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(X, k, n, mod);
        nmod_mat_init(X2, k, n, mod);
        nmod_mat_init(AX, m, n, mod);

        nmod_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1));
        nmod_mat_randtest(X2, state);

        nmod_mat_mul(B, A, X2);

        solved = nmod_mat_can_solve(X, A, B);

        nmod_mat_mul(AX, A, X);

        if (!solved || !nmod_mat_equal(B, AX))
            TEST_FUNCTION_FAIL(
                    "AX != B\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "X = %{nmod_mat}\n"
                    "AX = %{nmod_mat}\n",
                    A, B, X, AX);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(X);
        nmod_mat_clear(X2);
        nmod_mat_clear(AX);
    }

    TEST_FUNCTION_END(state);
}
