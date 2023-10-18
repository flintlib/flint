/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_entry, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t a;
        slong j, k;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        fmpz_mat_init(a, rows, cols);

        for (j = 0; j < rows; j++)
        {
            for (k = 0; k < cols; k++)
            {
                fmpz_set_ui(fmpz_mat_entry(a,j,k), 3*j + 7*k);
            }
        }

        for (j = 0; j < rows; j++)
        {
            for (k = 0; k < cols; k++)
            {
                if (fmpz_get_ui(fmpz_mat_entry(a,j,k)) != 3*j + 7*k)
                {
                    flint_printf("FAIL: get/set entry %wd,%wd\n", j, k);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mat_clear(a);
    }

    TEST_FUNCTION_END(state);
}
