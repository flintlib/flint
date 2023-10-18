/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "d_mat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(d_mat_equal, state)
{
    int i;

    /* check A != B if A, B have different dimensions
     * set A = B and check A == B
     * compare matrices with different entries */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        d_mat_t A, B, C, D, E;
        slong m, n, j;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        d_mat_init(A, m, n);
        d_mat_init(B, m, n);
        d_mat_init(C, m, n);
        d_mat_init(D, m + 1, n);
        d_mat_init(E, m, n + 1);

        if (d_mat_equal(A, D) || d_mat_equal(A, E))
        {
            flint_printf("FAIL: different dimensions should not be equal\n");
            fflush(stdout);
            flint_abort();
        }

        d_mat_randtest(A, state, 0, 0);
        d_mat_set(B, A);

        if (!d_mat_equal(A, B))
        {
            flint_printf("FAIL: copied matrices should be equal\n");
            fflush(stdout);
            flint_abort();
        }

        if (m && n)
        {
            j = n_randint(state, m * n);
            A->entries[j] += 1;

            if (d_mat_equal(A, B))
            {
                flint_printf("FAIL: modified matrices should not be equal\n");
                fflush(stdout);
                flint_abort();
            }
        }

        d_mat_clear(A);
        d_mat_clear(B);
        d_mat_clear(C);
        d_mat_clear(D);
        d_mat_clear(E);
    }

    TEST_FUNCTION_END(state);
}
