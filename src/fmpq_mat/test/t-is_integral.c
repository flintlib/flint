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
#include "fmpq.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_is_integral, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A;
        fmpz_mat_t B;
        slong n, m;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);

        fmpz_mat_randtest(B, state, 10);
        fmpq_mat_set_fmpz_mat(A, B);

        if (!fmpq_mat_is_integral(A))
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        if (m && n)
        {
            slong i, j;
            i = n_randint(state, m);
            j = n_randint(state, n);

            fmpq_set_si(fmpq_mat_entry(A, i, j), 1, 2);

            if (fmpq_mat_is_integral(A))
            {
                flint_printf("FAIL\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mat_clear(A);
        fmpz_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
