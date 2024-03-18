/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_is_zero, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        fmpz_mat_init(A, rows, cols);

        if (!fmpz_mat_is_zero(A))
        {
            flint_printf("FAIL!\n");
            fflush(stdout);
            flint_abort();
        }

        if (rows && cols)
        {
            fmpz_mat_randrank(A, state, FLINT_MIN(rows, cols), 100);
            if (fmpz_mat_is_zero(A))
            {
                flint_printf("FAIL!\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
