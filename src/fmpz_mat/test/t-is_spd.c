/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("is_spd....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: Gram matrices are positive definite iff full rank */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong n = n_randint(state, 10);
        slong rk;
        fmpz_mat_t A;
        int res;

        fmpz_mat_init(A, n, n);

        fmpz_mat_randtest(A, state, 1 + n_randint(state, 200));
        fmpz_mat_gram(A, A);
        rk = fmpz_mat_rank(A);
        res = fmpz_mat_is_spd(A);

        if ((rk < n && res) || (rk == n && !res))
        {
            flint_printf("FAIL\n");
            flint_printf("n = %wd, rk = %wd, res = %wd\n", n, rk, res);
            fmpz_mat_print_pretty(A);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(A);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
