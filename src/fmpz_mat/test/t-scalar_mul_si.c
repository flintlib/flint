/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_scalar_mul_si, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t C, B, A;
        slong rows, cols;
        slong c;

        rows = n_randint(state, 10);
        cols = n_randint(state, 10);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_mat_init(C, rows, cols);

        c = z_randtest(state);
        fmpz_mat_randtest(A, state, 100);

        fmpz_mat_scalar_mul_si(B, A, c);

        if (c == 0)
        {
            if (!fmpz_mat_is_zero(B))
            {
                flint_printf("FAIL!\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            fmpz_mat_scalar_divexact_si(C, B, c);

            if (!fmpz_mat_equal(C, A))
            {
                flint_printf("FAIL!\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
