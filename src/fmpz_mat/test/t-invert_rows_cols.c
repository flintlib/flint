/*
    Copyright (C) 2019 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_invert_rows_cols, state)
{
    slong n, rep;

    /* Rectangular inversion of rows and cols */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        slong i, j;
        fmpz_mat_t A, B;

        n = n_randint(state, 10);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, n);

        fmpz_mat_randtest(A, state, 1+n_randint(state, 3));

        fmpz_mat_set(B, A);
        fmpz_mat_invert_rows(A, NULL);
        fmpz_mat_invert_cols(A, NULL);

        for (i = 0; i < A->r; i++)
        {
            for (j =0; j < A->c; j++)
            {
                if (fmpz_cmp(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, A->r - i - 1, A->c - j - 1)) != 0)
                {
                    flint_printf("FAIL: B != A\n");
                    flint_printf("A:\n");
                    fmpz_mat_print_pretty(A);
                    flint_printf("B:\n");
                    fmpz_mat_print_pretty(B);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
