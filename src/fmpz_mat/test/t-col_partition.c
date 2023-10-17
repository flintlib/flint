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

TEST_FUNCTION_START(fmpz_mat_col_partition, state)
{
    int result = 0;
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B;
        slong m, n, j, k, p1 = 1, p2;
        slong * part;

        m = n_randint(state, 20) + 1;
        n = n_randint(state, 20) + 1;

        part = (slong *) flint_malloc(m*sizeof(slong));

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, m);

        /* set first row */
        for (j = 0; j < n; j++)
           fmpz_randtest(A->rows[0] + j, state, 100);
        /* ensure row is distinct */
        fmpz_set_ui(A->rows[0] + n_randint(state, n), 0);

        /* fill remaining rows */
        for (k = 1; k < m; k++)
        {
           if (n_randint(state, 2) == 0)
           {
              /* random row */
              for (j = 0; j < n; j++)
                 fmpz_randtest(A->rows[k] + j, state, 100);
              /* ensure row is distinct */
              fmpz_set_ui(A->rows[k] + n_randint(state, n), k);
              p1++;
           } else
           {
              /* same as last row */
              for (j = 0; j < n; j++)
                 fmpz_set(A->rows[k] + j, A->rows[k - 1] + j);
           }
        }

        /* swap random rows */
        for (k = 0; k < m; k++)
        {
           slong r1 = n_randint(state, m);
           slong r2 = n_randint(state, m);
           fmpz * t = A->rows[r1];
           A->rows[r1] = A->rows[r2];
           A->rows[r2] = t;
        }

        /* transpose so rows are now columns */
        fmpz_mat_transpose(B, A);

        /* check partition is correct */
        p2 = fmpz_mat_col_partition(part, B, 1);

        result = (p1 > n || p2 == p1);

        if (!result)
        {
            flint_printf("FAIL: col_partition failed\n");
            flint_printf("m = %ld, n = %ld, p1 = %ld, p2 = %ld\n", m, n, p1, p2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);

        flint_free(part);
    }

    TEST_FUNCTION_END(state);
}
