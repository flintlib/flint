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

TEST_FUNCTION_START(nmod_mat_solve_triu_classical, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, X, B, Y;
        mp_limb_t m;
        slong rows, cols;
        int unit;

        m = n_randtest_prime(state, 0);
        rows = n_randint(state, 100);
        cols = n_randint(state, 100);
        unit = n_randint(state, 2);

        nmod_mat_init(A, rows, rows, m);
        nmod_mat_init(B, rows, cols, m);
        nmod_mat_init(X, rows, cols, m);
        nmod_mat_init(Y, rows, cols, m);

        nmod_mat_randtriu(A, state, unit);
        nmod_mat_randtest(X, state);
        nmod_mat_mul(B, A, X);

        /* Check Y = A^(-1) * (A * X) = X */
        nmod_mat_solve_triu_classical(Y, A, B, unit);
        if (!nmod_mat_equal(Y, X))
        {
            flint_printf("FAIL!\n");
            flint_printf("A:\n");
            nmod_mat_print_pretty(A);
            flint_printf("X:\n");
            nmod_mat_print_pretty(X);
            flint_printf("B:\n");
            nmod_mat_print_pretty(B);
            flint_printf("Y:\n");
            nmod_mat_print_pretty(Y);
            fflush(stdout);
            flint_abort();
        }

        /* Check aliasing */
        nmod_mat_solve_triu_classical(B, A, B, unit);
        if (!nmod_mat_equal(B, X))
        {
            flint_printf("FAIL!\n");
            flint_printf("aliasing test failed");
            flint_printf("A:\n");
            nmod_mat_print_pretty(A);
            flint_printf("B:\n");
            nmod_mat_print_pretty(B);
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(X);
        nmod_mat_clear(Y);
    }

    TEST_FUNCTION_END(state);
}
