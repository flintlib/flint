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

TEST_FUNCTION_START(fmpz_mat_mul_multi_mod, state)
{
    fmpz_mat_t A, B, C, D;
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong m, n, k;

        m = n_randint(state, 100);
        n = n_randint(state, 100);
        k = n_randint(state, 100);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);

        flint_set_num_threads(1 + n_randint(state, 4));

        fmpz_mat_mul_classical(C, A, B);
        fmpz_mat_mul_multi_mod(D, A, B);

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

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong m, n, k;

        m = n_randint(state, 3);
        n = n_randint(state, 3);
        k = n_randint(state, 3);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        fmpz_mat_randtest(A, state, n_randint(state, 25000) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 25000) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);

        flint_set_num_threads(1 + n_randint(state, 4));

        fmpz_mat_mul_classical(C, A, B);
        fmpz_mat_mul_multi_mod(D, A, B);

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

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong n;

        n = n_randint(state, 100);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(C, n, n);
        fmpz_mat_init(D, n, n);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);

        flint_set_num_threads(1 + n_randint(state, 4));

        fmpz_mat_mul_classical(C, A, A);
        fmpz_mat_mul_multi_mod(D, A, A);

        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal (squaring)\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
