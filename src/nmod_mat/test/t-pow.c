/*
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_pow, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        mp_limb_t mod;
        slong m, j;
        ulong exp;

        mod = n_randtest_not_zero(state);
        m = n_randint(state, 20);
        exp = n_randint(state, 50);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, m, mod);
        nmod_mat_init(C, m, m, mod);
        nmod_mat_init(D, m, m, mod);
        nmod_mat_randtest(A, state);

        nmod_mat_pow(B, A, exp);
        nmod_mat_one(C);
        for(j = 1; j <= exp; j++)
        {
            nmod_mat_mul(D, C, A);
            nmod_mat_swap(D, C);
        }
        nmod_mat_pow(A, A, exp);

        if (!(nmod_mat_equal(C, B) && nmod_mat_equal(C, A)))
        {
            flint_printf("FAIL: results not equal\n");fflush(stdout);
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
