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

TEST_FUNCTION_START(nmod_mat_scalar_addmul_ui, state)
{
    slong m, n, mod, rep;
    mp_limb_t x;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C, D;

        m = n_randint(state, 40);
        n = n_randint(state, 40);
        mod = n_randtest_not_zero(state);
        x = n_randint(state, mod);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(C, m, n, mod);
        nmod_mat_init(D, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_scalar_mul(C, A, x);
        nmod_mat_add(C, B, C);
        nmod_mat_scalar_addmul_ui(D, B, A, x);

        if (!nmod_mat_equal(C, D))
        {
            flint_printf("FAIL\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
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
