/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2022 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_can_solve_dixon, state)
{
    int i;

    /* Solve random systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, X, AX;
        fmpz_mat_t M;
        fmpz_t den;
        int success;
        slong n, m, k, bits;

        n = n_randint(state, 10);
        m = n_randint(state, 10);
        k = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpz_init(den);

        fmpq_mat_init(A, n, k);
        fmpq_mat_init(B, n, m);
        fmpq_mat_init(X, k, m);
        fmpq_mat_init(AX, n, m);

        fmpz_mat_init(M, n, k);

        fmpz_mat_randrank(M, state, n_randint(state, FLINT_MIN(n, k) + 1), bits);
        if (i % 2)
            fmpz_mat_randops(M, state, n_randint(state, 2*n*k + 1));
        fmpz_randtest_not_zero(den, state, bits);
        fmpq_mat_set_fmpz_mat_div_fmpz(A, M, den);

        fmpq_mat_randtest(B, state, bits);

        success = fmpq_mat_can_solve_dixon(X, A, B);
        fmpq_mat_mul(AX, A, X);

        if (success && !fmpq_mat_equal(AX, B))
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

        fmpz_clear(den);

        fmpz_mat_clear(M);
        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(X);
        fmpq_mat_clear(AX);
    }

    /* Solve random soluble systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, X, AX;
        fmpz_mat_t M;
        fmpz_t den;
        int success;
        slong n, m, k, bits;

        n = n_randint(state, 10);
        m = n_randint(state, 10);
        k = n_randint(state, 10);

        bits = 1 + n_randint(state, 300);

        fmpz_init(den);

        fmpq_mat_init(A, n, k);
        fmpq_mat_init(B, n, m);
        fmpq_mat_init(X, k, m);
        fmpq_mat_init(AX, n, m);

        fmpz_mat_init(M, n, k);

        fmpz_mat_randrank(M, state, n_randint(state, FLINT_MIN(n, k) + 1), bits);
        if (i % 2)
            fmpz_mat_randops(M, state, n_randint(state, 2*m*n + 1));
        fmpz_randtest_not_zero(den, state, bits);
        fmpq_mat_set_fmpz_mat_div_fmpz(A, M, den);

        fmpq_mat_randtest(X, state, bits);

        fmpq_mat_mul(B, A, X);

        success = fmpq_mat_can_solve_dixon(X, A, B);
        fmpq_mat_mul(AX, A, X);

        if (!success || !fmpq_mat_equal(AX, B))
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

        fmpz_clear(den);

        fmpz_mat_clear(M);
        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(X);
        fmpq_mat_clear(AX);
    }

    TEST_FUNCTION_END(state);
}
