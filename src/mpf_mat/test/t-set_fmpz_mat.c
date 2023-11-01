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
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "mpf_mat.h"

TEST_FUNCTION_START(mpf_mat_set_fmpz_mat, state)
{
    int i;

    /* set entries of an fmpz_mat, convert to mpf_mat and then check that
       the entries remain same */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A;
        mpf_mat_t B;
        slong j, k;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        fmpz_mat_init(A, rows, cols);
        mpf_mat_init(B, rows, cols, mpf_get_default_prec());

        for (j = 0; j < rows; j++)
        {
            for (k = 0; k < cols; k++)
            {
                fmpz_set_ui(fmpz_mat_entry(A, j, k), 3 * j + 7 * k);
            }
        }

        mpf_mat_set_fmpz_mat(B, A);

        for (j = 0; j < rows; j++)
        {
            for (k = 0; k < cols; k++)
            {
                if (flint_mpf_cmp_ui(mpf_mat_entry(B, j, k), 3 * j + 7 * k) != 0)
                {
                    flint_printf("FAIL: j = %wd, k = %wd\n", j, k);
                    fmpz_mat_print_pretty(A);
                    mpf_mat_print(B);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mat_clear(A);
        mpf_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
