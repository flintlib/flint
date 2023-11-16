/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_exp_taylor_sum, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_mat_t A, S1, S2, S3;
        fmpz_t f;
        slong n, N, prec1, prec2;

        n = n_randint(state, 5);
        N = n_randint(state, 40);
        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);

        arb_mat_init(A, n, n);
        arb_mat_init(S1, n, n);
        arb_mat_init(S2, n, n);
        arb_mat_init(S3, n, n);
        fmpz_init(f);

        arb_mat_randtest(A, state, prec1, 10);
        arb_mat_randtest(S1, state, prec1, 10);
        arb_mat_randtest(S2, state, prec1, 10);

        arb_mat_exp_taylor_sum(S1, A, N, prec1);

        arb_mat_exp_taylor_sum(S2, A, N + 1, prec2);

        arb_mat_pow_ui(S3, A, N, prec2);
        fmpz_fac_ui(f, N);
        arb_mat_scalar_div_fmpz(S3, S3, f, prec2);
        arb_mat_add(S3, S3, S1, prec2);

        if (!arb_mat_overlaps(S2, S3))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = %wd, N = %wd\n", n, N);
            flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("S1 = \n"); arb_mat_printd(S1, 15); flint_printf("\n\n");
            flint_printf("S2 = \n"); arb_mat_printd(S2, 15); flint_printf("\n\n");
            flint_printf("S3 = \n"); arb_mat_printd(S3, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_mat_exp_taylor_sum(A, A, N, prec1);

        if (!arb_mat_overlaps(A, S1))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(S1);
        arb_mat_clear(S2);
        arb_mat_clear(S3);
        fmpz_clear(f);
    }

    TEST_FUNCTION_END(state);
}
