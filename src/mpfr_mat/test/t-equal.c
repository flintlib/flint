/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpfr_mat.h"

TEST_FUNCTION_START(mpfr_mat_equal, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpfr_mat_t A, B, D, E;
        slong m, n, j;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mpfr_mat_init(A, m, n, 200);
        mpfr_mat_init(B, m, n, 200);
        mpfr_mat_init(D, m + 1, n, 200);
        mpfr_mat_init(E, m, n + 1, 200);

        if (mpfr_mat_equal(A, D) || mpfr_mat_equal(A, E))
            TEST_FUNCTION_FAIL("different dimensions should not be equal\n");

        mpfr_mat_randtest(A, state);
        mpfr_mat_set(B, A);

        if (!mpfr_mat_equal(A, B))
            TEST_FUNCTION_FAIL("copied matrices should be equal\n");

        if (m && n)
        {
            j = n_randint(state, m * n);
            mpfr_add_ui(A->entries + j, A->entries + j, 1, MPFR_RNDN);

            if (mpfr_mat_equal(A, B))
                TEST_FUNCTION_FAIL("modified matrices should not be equal\n");
        }

        mpfr_mat_clear(A);
        mpfr_mat_clear(B);
        mpfr_mat_clear(D);
        mpfr_mat_clear(E);
    }

    TEST_FUNCTION_END(state);
}
