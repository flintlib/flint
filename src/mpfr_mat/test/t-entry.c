/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpfr_mat.h"

TEST_FUNCTION_START(mpfr_mat_entry, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpfr_mat_t A;
        slong j, k;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        mpfr_mat_init(A, rows, cols, 200);

        for (j = 0; j < rows; j++)
        {
            for (k = 0; k < cols; k++)
            {
                mpfr_set_si(mpfr_mat_entry(A, j, k), 3 * j + 7 * k, MPFR_RNDN);
            }
        }

        for (j = 0; j < rows; j++)
        {
            for (k = 0; k < cols; k++)
            {
                if (mpfr_cmp_si(mpfr_mat_entry(A, j, k), 3 * j + 7 * k) != 0)
                {
                    flint_printf("FAIL: get/set entry %wd, %wd\n", j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        mpfr_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
