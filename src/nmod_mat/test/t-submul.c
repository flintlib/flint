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

TEST_FUNCTION_START(nmod_mat_submul, state)
{
    slong i;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D, T, E;
        mp_limb_t mod = n_randtest_not_zero(state);

        slong m, k, n;

        m = n_randint(state, 100);
        k = n_randint(state, 100);
        n = n_randint(state, 100);

        /* Force Strassen test */
        if (i < 5)
        {
            m += 300;
            k += 300;
            n += 300;
        }

        nmod_mat_init(A, m, k, mod);
        nmod_mat_init(B, k, n, mod);
        nmod_mat_init(C, m, n, mod);
        nmod_mat_init(D, m, n, mod);
        nmod_mat_init(T, m, n, mod);
        nmod_mat_init(E, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);
        nmod_mat_randtest(C, state);

        nmod_mat_submul(D, C, A, B);

        nmod_mat_mul(T, A, B);
        nmod_mat_sub(E, C, T);

        if (!nmod_mat_equal(D, E))
        {
            flint_printf("FAIL: results not equal\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
            nmod_mat_print_pretty(E);
            fflush(stdout);
            flint_abort();
        }

        /* Check aliasing */
        nmod_mat_submul(C, C, A, B);

        if (!nmod_mat_equal(C, E))
        {
            flint_printf("FAIL: results not equal (aliasing)\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
            nmod_mat_print_pretty(E);
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
        nmod_mat_clear(E);
        nmod_mat_clear(T);
    }

    TEST_FUNCTION_END(state);
}
