/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_scalar_mul_2exp, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t C, B, A;
        slong rows, cols;
        ulong exp;

        rows = n_randint(state, 10);
        cols = n_randint(state, 10);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_mat_init(C, rows, cols);

        exp = n_randint(state, 200);
        fmpz_mat_randtest(A, state, 100);

        fmpz_mat_scalar_mul_2exp(B, A, exp);
        fmpz_mat_scalar_tdiv_q_2exp(C, B, exp);

        if (!fmpz_mat_equal(C, A))
        {
            flint_printf("FAIL!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    /* test aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t C, B, A;
        slong rows, cols;
        ulong exp;

        rows = n_randint(state, 10);
        cols = n_randint(state, 10);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_mat_init(C, rows, cols);

        exp = n_randint(state, 200);
        fmpz_mat_randtest(A, state, 100);

        fmpz_mat_scalar_mul_2exp(B, A, exp);
        fmpz_mat_scalar_tdiv_q_2exp(C, B, exp);
        fmpz_mat_scalar_mul_2exp(A, A, exp);
        fmpz_mat_scalar_tdiv_q_2exp(A, A, exp);

        if (!fmpz_mat_equal(C, A))
        {
            flint_printf("FAIL!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
