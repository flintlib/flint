/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_gram, state)
{
	fmpz_mat_t A, B, C, D;
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong m, n;

        m = n_randint(state, 50);
        n = n_randint(state, 50);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, m);
        fmpz_mat_init(C, m, m);
        fmpz_mat_init(D, m, m);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);

        fmpz_mat_transpose(B, A);
        fmpz_mat_mul(C, A, B);
        fmpz_mat_gram(D, A);

        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
