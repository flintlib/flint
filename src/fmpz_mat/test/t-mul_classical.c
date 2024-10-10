/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_mul_classical, state)
{
    fmpz_mat_t A, B, C, D;
    slong ix;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        slong m, n, k;
        slong i, j, h;

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        k = n_randint(state, 50);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);

        fmpz_mat_mul_classical(C, A, B);

        for (i = 0; i < D->r; i++)
            for (j = 0; j < D->c; j++)
                for (h = 0; h < B->r; h++)
                    fmpz_addmul(fmpz_mat_entry(D, i, j),
                        fmpz_mat_entry(A, i, h),
                        fmpz_mat_entry(B, h, j));

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
