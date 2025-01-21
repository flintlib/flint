/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_init_clear, state)
{
    slong m, n, mod, rep;

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A;

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);

        if (!nmod_mat_is_zero(A))
            TEST_FUNCTION_FAIL("entries not zero\n");

        if (A->mod.n != mod)
            TEST_FUNCTION_FAIL("bad modulus\n");

        nmod_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
