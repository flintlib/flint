/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_transpose, state)
{
    slong m, n, mod, mod2, rep;

    /* Rectangular transpose, same modulus */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, m, mod);
        nmod_mat_init(C, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_transpose(B, A);
        nmod_mat_transpose(C, B);

        if (!nmod_mat_equal(C, A))
        {
            flint_printf("FAIL: C != A\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
    }

    /* Rectangular transpose, different modulus */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, AT, B, BT, AT2;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mod = n_randtest_not_zero(state);
        mod2 = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(AT, n, m, mod);
        nmod_mat_init(B, m, n, mod2);
        nmod_mat_init(BT, n, m, mod2);
        nmod_mat_init(AT2, n, m, mod2);

        nmod_mat_randtest(A, state);
        nmod_mat_set(B, A);

        nmod_mat_transpose(AT, A);
        nmod_mat_transpose(BT, B);

        nmod_mat_set(AT2, AT);

        if (!nmod_mat_equal(BT, AT2))
        {
            flint_printf("FAIL: AT != BT\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(AT);
        nmod_mat_clear(AT2);
        nmod_mat_clear(B);
        nmod_mat_clear(BT);
    }

    /* Self-transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B;

        m = n_randint(state, 20);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, m, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_set(B, A);

        nmod_mat_transpose(B, B);
        nmod_mat_transpose(B, B);

        if (!nmod_mat_equal(B, A))
        {
            flint_printf("FAIL: B != A\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
