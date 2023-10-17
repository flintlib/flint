/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_mul_small, state)
{
    fmpz_mat_t A, B, C, D;
    slong i;
    slong max_threads = 6;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong m, k, n;

        if (flint_get_num_threads() >= max_threads - 1)
        {
            m = n_randint(state, 200);
            n = n_randint(state, 200);
            k = n_randint(state, 200);
        }
        else
        {
            m = n_randint(state, 50);
            n = n_randint(state, 50);
            k = n_randint(state, 50);
        }

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, k, n);
        fmpz_mat_init(C, m, n);
        fmpz_mat_init(D, m, n);

        fmpz_mat_randtest(A, state, n_randint(state, SMALL_FMPZ_BITCOUNT_MAX) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, SMALL_FMPZ_BITCOUNT_MAX) + 1);
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(D, state, n_randint(state, 200) + 1);

        _fmpz_mat_mul_small(C, A, B);
        fmpz_mat_mul_classical_inline(D, A, B);

        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n\n");
            fmpz_mat_print(A); flint_printf("\n\n");
            fmpz_mat_print(B); flint_printf("\n\n");
            fmpz_mat_print(C); flint_printf("\n\n");
            fmpz_mat_print(D); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
