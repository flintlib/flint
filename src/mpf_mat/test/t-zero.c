/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "mpf_mat.h"

TEST_FUNCTION_START(mpf_mat_zero, state)
{
    slong m, n, i, j, rep;

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        mpf_mat_t A;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mpf_mat_init(A, m, n, 200);

        mpf_mat_randtest(A, state, 200);
        mpf_mat_zero(A);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (flint_mpf_cmp_ui(mpf_mat_entry(A, i, j), 0) != 0)
                {
                    flint_printf("FAIL: nonzero entry\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        mpf_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
