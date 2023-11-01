/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_scalar_smod, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t B, A;
        slong j, k, rows, cols;
        fmpz_t P, c;

        rows = n_randint(state, 10);
        cols = n_randint(state, 10);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_init(P);
        fmpz_init(c);

        do {
           fmpz_randtest(P, state, 100);
           fmpz_abs(P, P);
        } while (fmpz_is_zero(P));

        fmpz_mat_randtest(A, state, 100);

        fmpz_mat_scalar_smod(B, A, P);

        for (j = 0; j < A->r; j++)
        {
           for (k = 0; k < A->c; k++)
           {
              fmpz_smod(c, A->rows[j] + k, P);

              if (!fmpz_equal(c, B->rows[j] + k))
              {
                 flint_printf("FAIL!\n");
                 flint_printf("%wd, %wd\n", j, k);
                 fmpz_print(c); flint_printf("\n");
                 fmpz_print(A->rows[j] + k); flint_printf("\n");
                 fmpz_print(B->rows[j] + k); flint_printf("\n");
                 fmpz_print(P); flint_printf("\n");
                 fflush(stdout);
                 flint_abort();
              }
           }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_clear(P);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
