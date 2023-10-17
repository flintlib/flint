/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_mul_blas, state)
{
    slong i;
    slong max_num_threads = 5;

    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        slong m, n, k;
        fmpz_mat_t A, B, C, D;

        m = 1 + n_randint(state, 50);
        n = 1 + n_randint(state, 50);
        k = 1 + n_randint(state, 50);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        if (n_randint(state, 2))
            fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        else
            fmpz_mat_randtest_unsigned(A, state, n_randint(state, 200) + 1);

        if (n_randint(state, 2))
            fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);
        else
            fmpz_mat_randtest_unsigned(B, state, n_randint(state, 200) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);

        fmpz_mat_mul_classical_inline(C, A, B);
        if (fmpz_mat_mul_blas(D, A, B))
        {
            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal\n");
                fflush(stdout);
                flint_abort();
            }
        }
#if FLINT_USES_BLAS && FLINT_BITS == 64
        else
        {
            flint_printf("FAIL: blas should have worked\n");
            fflush(stdout);
            flint_abort();
        }
#endif

        flint_set_num_threads(1 + n_randint(state, max_num_threads));

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        slong m, n, k;
        fmpz_mat_t A, B, C, D;

        m = 1 + n_randint(state, 3);
        n = 1 + n_randint(state, 3);
        k = 1 + n_randint(state, 3);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        fmpz_mat_randtest(A, state, n_randint(state, 20000) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 20000) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(D, state, n_randint(state, 200) + 1);

        fmpz_mat_mul_classical_inline(C, A, B);
        if (fmpz_mat_mul_blas(D, A, B))
        {
            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal\n");
                fflush(stdout);
                flint_abort();
            }
        }
#if FLINT_USES_BLAS && FLINT_BITS == 64
        else
        {
            flint_printf("FAIL: blas should have worked\n");
            fflush(stdout);
            flint_abort();
        }
#endif

        flint_set_num_threads(1 + n_randint(state, max_num_threads));

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
