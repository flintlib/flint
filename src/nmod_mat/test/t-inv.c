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

TEST_FUNCTION_START(nmod_mat_inv, state)
{
    nmod_mat_t A, B, C, I;
    slong i, j, m, r;
    mp_limb_t mod;
    int result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, m, mod);
        nmod_mat_init(C, m, m, mod);
        nmod_mat_init(I, m, m, mod);

        for (j = 0; j < m; j++)
            I->rows[j][j] = UWORD(1);

        /* Verify that A * A^-1 = I for random matrices */

        nmod_mat_randrank(A, state, m);
        /* Dense or sparse? */
        if (n_randint(state, 2))
            nmod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        result = nmod_mat_inv(B, A);
        nmod_mat_mul(C, A, B);

        if (!nmod_mat_equal(C, I) || !result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A * A^-1 != I!\n");
            flint_printf("A:\n");
            nmod_mat_print_pretty(A);
            flint_printf("A^-1:\n");
            nmod_mat_print_pretty(B);
            flint_printf("A * A^-1:\n");
            nmod_mat_print_pretty(C);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* Test aliasing */
        nmod_mat_set(C, A);
        nmod_mat_inv(A, A);
        nmod_mat_mul(B, A, C);

        if (!nmod_mat_equal(B, I))
        {
            flint_printf("FAIL:\n");
            flint_printf("aliasing failed!\n");
            nmod_mat_print_pretty(C);
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(I);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        mod = n_randtest_prime(state, 0);
        r = n_randint(state, m);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, m, mod);

        nmod_mat_randrank(A, state, r);

        /* Dense */
        if (n_randint(state, 2))
            nmod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        result = nmod_mat_inv(B, A);

        if (result)
        {
            flint_printf("FAIL:\n");
            flint_printf("singular matrix reported as invertible\n");
            fflush(stdout);
            flint_abort();
        }

        /* Aliasing */
        result = nmod_mat_inv(A, A);
        if (result)
        {
            flint_printf("FAIL:\n");
            flint_printf("singular matrix reported as invertible\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
