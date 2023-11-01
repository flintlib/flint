/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

#define FMPZ_MAT_CHOL_EPS (1.0E-9)

TEST_FUNCTION_START(fmpz_mat_chol_d, state)
{
    int i;

    /* check RR^T = A */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A;
        d_mat_t R, Rt, Atmp, Btmp;

        slong m;

        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);

        d_mat_init(R, m, m);
        d_mat_init(Rt, m, m);
        d_mat_init(Atmp, m, m);
        d_mat_init(Btmp, m, m);

        fmpz_mat_randtest(A, state, 10);
        fmpz_mat_gram(A, A);
        fmpz_mat_get_d_mat(Atmp, A);
        d_mat_zero(R);

        fmpz_mat_chol_d(R, A);
        d_mat_transpose(Rt, R);

        d_mat_mul_classical(Btmp, R, Rt);

        if (!d_mat_approx_equal(Atmp, Btmp, FMPZ_MAT_CHOL_EPS))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpz_mat_print_pretty(A);
            flint_printf("R:\n");
            d_mat_print(R);
            flint_printf("R^T:\n");
            d_mat_print(Rt);
            flint_printf("Btmp:\n");
            d_mat_print(Btmp);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);

        d_mat_clear(R);
        d_mat_clear(Rt);
        d_mat_clear(Atmp);
        d_mat_clear(Btmp);
    }

    TEST_FUNCTION_END(state);
}
