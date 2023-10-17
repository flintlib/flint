/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_add_sub, state)
{
    slong m, n, rep;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mat_t A;
        fmpz_mat_t B;
        fmpz_mat_t C;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(C, m, n);

        fmpz_mat_randtest(A, state, 100);
        fmpz_mat_randtest(B, state, 100);

        fmpz_mat_neg(C, A);
        fmpz_mat_add(A, A, B);
        fmpz_mat_sub(A, A, B);
        fmpz_mat_neg(A, A);

        if (!fmpz_mat_equal(A, C))
        {
            flint_printf("FAIL: matrices not equal!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
