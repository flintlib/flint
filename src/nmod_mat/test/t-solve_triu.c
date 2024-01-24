/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_solve_triu, state)
{
    slong i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, X, B, Y;
        mp_limb_t m;
        slong rows, cols;
        int unit;

        m = n_randtest_prime(state, 0);
        rows = n_randint(state, 200);
        cols = n_randint(state, 200);
        unit = n_randint(state, 2);

        nmod_mat_init(A, rows, rows, m);
        nmod_mat_init(B, rows, cols, m);
        nmod_mat_init(X, rows, cols, m);
        nmod_mat_init(Y, rows, cols, m);

        nmod_mat_randtriu(A, state, unit);
        nmod_mat_randtest(X, state);
        nmod_mat_mul(B, A, X);

        /* Check Y = A^(-1) * (A * X) = X */
        nmod_mat_solve_triu(Y, A, B, unit);
        if (!nmod_mat_equal(Y, X))
            TEST_FUNCTION_FAIL(
                    "A = %{nmod_mat}\n"
                    "X = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "Y = %{nmod_mat}\n",
                    A, X, B, Y);

        /* Check aliasing */
        nmod_mat_solve_triu(B, A, B, unit);
        if (!nmod_mat_equal(B, X))
            TEST_FUNCTION_FAIL(
                    "Aliasing test failed\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n",
                    A, B);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(X);
        nmod_mat_clear(Y);
    }

    TEST_FUNCTION_END(state);
}
