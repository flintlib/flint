/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_solve_vec, state)
{
    nmod_mat_t A, x, b, Ax;
    slong i, m, r;
    int solved;
    mp_limb_t mod;

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(b, m, 1, mod);
        nmod_mat_init(x, m, 1, mod);
        nmod_mat_init(Ax, m, 1, mod);

        nmod_mat_randrank(A, state, m);
        nmod_mat_randtest(b, state);

        /* Dense */
        if (n_randint(state, 2))
            nmod_mat_randops(A, state, 1+n_randint(state, 1+m*m));

        solved = nmod_mat_solve_vec(x->entries, A, b->entries);
        nmod_mat_mul(Ax, A, x);

        if (!nmod_mat_equal(Ax, b) || !solved)
            TEST_FUNCTION_FAIL(
                    "Ax != b\n"
                    "A = %{nmod_mat}\n"
                    "b = %{nmod_mat}\n"
                    "x = %{nmod_mat}\n"
                    "Ax = %{nmod_mat}\n",
                    A, b, x, Ax);

        nmod_mat_clear(A);
        nmod_mat_clear(b);
        nmod_mat_clear(x);
        nmod_mat_clear(Ax);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        r = n_randint(state, m);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(b, m, 1, mod);
        nmod_mat_init(x, m, 1, mod);
        nmod_mat_init(Ax, m, 1, mod);

        nmod_mat_randrank(A, state, r);
        nmod_mat_randtest(b, state);

        /* Dense */
        if (n_randint(state, 2))
            nmod_mat_randops(A, state, 1+n_randint(state, 1+m*m));

        solved = nmod_mat_solve_vec(x->entries, A, b->entries);

        if (solved)
            TEST_FUNCTION_FAIL("singular system was 'solved'\n");

        nmod_mat_clear(A);
        nmod_mat_clear(b);
        nmod_mat_clear(x);
        nmod_mat_clear(Ax);
    }

    TEST_FUNCTION_END(state);
}
