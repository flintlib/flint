/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_solve, state)
{
    nmod_mat_t A, X, B, AX;
    slong i, m, n, r;
    mp_limb_t mod;
    int solved;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(X, m, n, mod);
        nmod_mat_init(AX, m, n, mod);

        nmod_mat_randrank(A, state, m);
        nmod_mat_randtest(B, state);

        /* Dense */
        if (n_randint(state, 2))
            nmod_mat_randops(A, state, 1+n_randint(state, 1+m*m));

        solved = nmod_mat_solve(X, A, B);

        nmod_mat_mul(AX, A, X);

        if (!nmod_mat_equal(AX, B) || !solved)
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

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);
        r = n_randint(state, m);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(X, m, n, mod);
        nmod_mat_init(AX, m, n, mod);

        nmod_mat_randrank(A, state, r);
        nmod_mat_randtest(B, state);

        /* Dense */
        if (n_randint(state, 2))
            nmod_mat_randops(A, state, 1+n_randint(state, 1+m*m));

        solved = nmod_mat_solve(X, A, B);

        if (solved)
            TEST_FUNCTION_FAIL(
                    "singular system was 'solved'\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "X = %{nmod_mat}\n",
                    A, B, X);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(X);
        nmod_mat_clear(AX);
    }

    TEST_FUNCTION_END(state);
}
