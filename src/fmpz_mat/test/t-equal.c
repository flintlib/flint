/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_equal, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C, D, E;
        slong m, n, j;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(C, m, n);
        fmpz_mat_init(D, m+1, n);
        fmpz_mat_init(E, m, n+1);

        if (fmpz_mat_equal(A, D) || fmpz_mat_equal(A, E))
        {
            flint_printf("FAIL: different dimensions should not be equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_randtest(A, state, 1 + n_randint(state, 100));
        fmpz_mat_set(B, A);

        if (!fmpz_mat_equal(A, B))
        {
            flint_printf("FAIL: copied matrices should be equal\n");
            fflush(stdout);
            flint_abort();
        }

        if (m && n)
        {
            j = n_randint(state, m*n);
            fmpz_add_ui(A->entries + j, A->entries + j, 1);

            if (fmpz_mat_equal(A, B))
            {
                flint_printf("FAIL: modified matrices should not be equal\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(E);
    }

    TEST_FUNCTION_END(state);
}
