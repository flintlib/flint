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

TEST_FUNCTION_START(nmod_mat_init_clear, state)
{
    slong m, n, mod, i, j, rep;

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A;

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (A->rows[i][j] != UWORD(0))
                {
                    flint_printf("FAIL: entries not zero!\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        if (A->mod.n != mod)
        {
            flint_printf("FAIL: bad modulus\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
