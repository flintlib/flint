/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_pow, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        slong i, n;
        ulong e;

        n = n_randint(state, 10);
        e = n_randint(state, 20);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, n);
        fmpz_mat_init(C, n, n);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        fmpz_mat_pow(B, A, e);

        fmpz_mat_one(C);
        for (i = 0; i < e; i++)
            fmpz_mat_mul(C, C, A);

        if (!fmpz_mat_equal(C, B))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_pow(A, A, e);

        if (!fmpz_mat_equal(A, B))
        {
            flint_printf("FAIL: aliasing failed\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
