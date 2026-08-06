/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_is_diagonal, state)
{
    slong iter;

    /* Zero matrix is diagonal */
    {
        fmpz_mat_t A;
        fmpz_mat_init(A, 5, 5);
        fmpz_mat_zero(A);
        if (!fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: zero matrix\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }

    /* Identity matrix is diagonal */
    {
        fmpz_mat_t A;
        fmpz_mat_init(A, 4, 4);
        fmpz_mat_one(A);
        if (!fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: identity matrix\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }

    /* 0x0 matrix is diagonal */
    {
        fmpz_mat_t A;
        fmpz_mat_init(A, 0, 0);
        if (!fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: 0x0 matrix\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }

    /* 0x5 matrix is diagonal */
    {
        fmpz_mat_t A;
        fmpz_mat_init(A, 0, 5);
        if (!fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: 0x5 matrix\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }

    /* 3x0 matrix is diagonal */
    {
        fmpz_mat_t A;
        fmpz_mat_init(A, 3, 0);
        if (!fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: 3x0 matrix\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }

    /* Non-square diagonal matrix */
    {
        fmpz_mat_t A;
        fmpz_mat_init(A, 3, 5);
        fmpz_mat_zero(A);
        fmpz_set_si(fmpz_mat_entry(A, 0, 0), 7);
        fmpz_set_si(fmpz_mat_entry(A, 1, 1), -3);
        fmpz_set_si(fmpz_mat_entry(A, 2, 2), 11);
        if (!fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: 3x5 diagonal matrix\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }

    /* Non-diagonal matrix */
    {
        fmpz_mat_t A;
        fmpz_mat_init(A, 3, 3);
        fmpz_mat_zero(A);
        fmpz_set_si(fmpz_mat_entry(A, 0, 0), 1);
        fmpz_set_si(fmpz_mat_entry(A, 0, 1), 2);
        if (fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: should not be diagonal\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }

    /* 1x1 matrix is always diagonal */
    {
        fmpz_mat_t A;
        fmpz_mat_init(A, 1, 1);
        fmpz_set_si(fmpz_mat_entry(A, 0, 0), -42);
        if (!fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: 1x1 matrix\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }

    /* Randomized: diagonal matrices should return 1 */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A;
        slong m, n, d, i;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        d = FLINT_MIN(m, n);

        fmpz_mat_init(A, m, n);
        fmpz_mat_zero(A);
        for (i = 0; i < d; i++)
            fmpz_randtest(fmpz_mat_entry(A, i, i), state, 100);

        if (!fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: random diagonal m=%wd n=%wd\n", m, n);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
    }

    /* Randomized: non-diagonal matrices should return 0 */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A;
        slong m, n, i, j;

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        if (m == 1 && n == 1)
            continue;

        fmpz_mat_init(A, m, n);
        fmpz_mat_zero(A);

        /* Place a nonzero entry off-diagonal */
        do
        {
            i = n_randint(state, m);
            j = n_randint(state, n);
        } while (i == j);

        fmpz_set_si(fmpz_mat_entry(A, i, j),
            1 + (slong) n_randint(state, 100));

        if (fmpz_mat_is_diagonal(A))
        {
            flint_printf("FAIL: should not be diagonal, "
                "nonzero at (%wd,%wd)\n", i, j);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
