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

TEST_FUNCTION_START(nmod_mat_invert_rows_cols, state)
{
    slong m, n, mod, rep;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B;
        slong i, j;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);

        nmod_mat_randtest(A, state);

        nmod_mat_set(B, A);

        nmod_mat_invert_rows(A, NULL);
        nmod_mat_invert_cols(A, NULL);

        for (i = 0; i < A->r; i++)
        {
            for (j =0; j < A->c; j++)
            {
                if (nmod_mat_entry(B, i, j) != nmod_mat_entry(A, A->r - i - 1, A->c - j - 1))
                {
                    flint_printf("FAIL: B != A\n");
                    flint_printf("A:\n");
                    nmod_mat_print_pretty(A);
                    flint_printf("B:\n");
                    nmod_mat_print_pretty(B);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
