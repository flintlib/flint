/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_is_spd, state)
{
    slong iter;

    /* Test:
       - fails for non-square matrices
       - fails for non-symmetric matrices
       - Gram matrices are positive definite iff full rank */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong n = n_randint(state, 10);
        slong rk;
        fmpz_mat_t A, B, C;
        slong j, k;
        int res;

        fmpz_mat_init(A, n, n + 1 + n_randint(state, 10));
        fmpz_mat_init(B, n, n);
        fmpz_mat_init(C, n, n);

        fmpz_mat_one(A);
        if (fmpz_mat_is_spd(A))
        {
            flint_printf("FAIL (square)\n");
            flint_abort();
        }

        fmpz_mat_randtest(B, state, 1 + n_randint(state, 200));
        fmpz_mat_gram(B, B);

        if (n > 1)
        {
            fmpz_mat_set(C, B);
            j = n_randint(state, n - 1);
            k = j + 1 + n_randint(state, n - j - 1);
            fmpz_add_si(fmpz_mat_entry(C, j, k), fmpz_mat_entry(C, j, k), 1);
            if (fmpz_mat_is_spd(C))
            {
                flint_printf("FAIL (symmetric)\n");
                flint_abort();
            }
        }

        rk = fmpz_mat_rank(B);
        res = fmpz_mat_is_spd(B);

        if ((rk < n && res) || (rk == n && !res))
        {
            flint_printf("FAIL\n");
            flint_printf("n = %wd, rk = %wd, res = %wd\n", n, rk, res);
            fmpz_mat_print_pretty(B);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
