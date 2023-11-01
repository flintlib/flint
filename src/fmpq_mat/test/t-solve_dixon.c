/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_solve_dixon, state)
{
    int i;

    /* Solve nonsingular systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, X, AX;
        fmpq_t d;
        int success;
        slong n, m, bits;

        n = n_randint(state, 10);
        m = n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, m);
        fmpq_mat_init(X, n, m);
        fmpq_mat_init(AX, n, m);

        fmpq_init(d);
        /* XXX: replace with a randtest function */
        do {
            fmpq_mat_randtest(A, state, bits);
            fmpq_mat_det(d, A);
        } while (fmpq_is_zero(d));
        fmpq_clear(d);

        fmpq_mat_randtest(B, state, bits);

        success = fmpq_mat_solve_dixon(X, A, B);
        fmpq_mat_mul(AX, A, X);

        if (!fmpq_mat_equal(AX, B) || !success)
        {
            flint_printf("FAIL!\n");
            flint_printf("success: %d\n", success);
            flint_printf("A:\n");
            fmpq_mat_print(A);
            flint_printf("B:\n");
            fmpq_mat_print(B);
            flint_printf("X:\n");
            fmpq_mat_print(X);
            flint_printf("AX:\n");
            fmpq_mat_print(AX);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(X);
        fmpq_mat_clear(AX);
    }

    /* Check singular systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, X;
        fmpz_mat_t M;
        fmpz_t den;
        slong n, m, bits;
        int success;

        n = 1 + n_randint(state, 10);
        m = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_init(den);
        fmpz_mat_init(M, n, n);
        fmpz_mat_randrank(M, state, n_randint(state, n), bits);
        if (i % 2)
            fmpz_mat_randops(M, state, n_randint(state, 2*m*n + 1));
        fmpz_randtest_not_zero(den, state, bits);
        fmpq_mat_init(A, n, n);
        fmpq_mat_set_fmpz_mat_div_fmpz(A, M, den);

        fmpq_mat_init(B, n, m);
        fmpq_mat_randtest(B, state, bits);
        fmpq_mat_init(X, n, m);

        success = fmpq_mat_solve_dixon(X, A, B);

        if (success != 0)
        {
            flint_printf("FAIL!\n");
            flint_printf("Expected success = 0\n");
            fmpq_mat_print(A);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(X);
        fmpz_mat_clear(M);
        fmpz_clear(den);
    }

    TEST_FUNCTION_END(state);
}
