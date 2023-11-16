/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_spd_get_fmpz_mat, state)
{
    slong iter;

    /* Test: construct input from an integral matrix */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong n = n_randint(state, 10);
        slong prec = 200;
        slong mag_exp = -prec - 10;
        slong rk = 0;
        fmpz_mat_t N, T;
        arb_mat_t A;
        mag_t err;
        int res;

        fmpz_mat_init(N, n, n);
        fmpz_mat_init(T, n, n);
        arb_mat_init(A, n, n);
        mag_init(err);

        while (rk < n)
        {
            fmpz_mat_randtest(N, state, 1 + n_randint(state, 50));
            fmpz_mat_gram(N, N);
            rk = fmpz_mat_rank(N);
        }

        arb_mat_set_fmpz_mat(A, N);
        arb_mat_scalar_mul_2exp_si(A, A, -prec);
        mag_one(err);
        mag_mul_2exp_si(err, err, mag_exp);
        arb_mat_add_error_mag(A, err);

        res = arb_mat_spd_get_fmpz_mat(T, A, prec);

        if (!res || !fmpz_mat_equal(T, N))
        {
            flint_printf("FAIL (res = %wd)\n", res);
            fmpz_mat_print_pretty(T);
            flint_printf("\n");
            fmpz_mat_print_pretty(N);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(N);
        fmpz_mat_clear(T);
        arb_mat_clear(A);
        mag_clear(err);
    }

    TEST_FUNCTION_END(state);
}
