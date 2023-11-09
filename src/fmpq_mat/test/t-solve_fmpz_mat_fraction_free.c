/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

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

TEST_FUNCTION_START(fmpq_mat_solve_fmpz_mat_fraction_free, state)
{
    int i;

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t X, AX;
        fmpz_mat_t A, B, AX_Z;
        fmpz_t d;
        int success;

        slong n, m, bits;

        n = n_randint(state, 40);
        m = n_randint(state, 40);
        bits = 1 + n_randint(state, 100);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, m);
        fmpz_mat_init(AX_Z, n, m);
        fmpq_mat_init(X, n, m);
        fmpq_mat_init(AX, n, m);

        fmpz_init(d);
        /* XXX: replace with a randtest function */
        do {
            fmpz_mat_randtest(A, state, bits);
            fmpz_mat_det(d, A);
        } while (fmpz_is_zero(d));
        fmpz_clear(d);

        fmpz_mat_randtest(B, state, bits);

        success = fmpq_mat_solve_fmpz_mat_fraction_free(X, A, B);
        fmpq_mat_mul_r_fmpz_mat(AX, A, X);

        if (!success || !fmpq_mat_get_fmpz_mat(AX_Z, AX)
                || !fmpz_mat_equal(AX_Z, B))
        {
            flint_printf("FAIL!\n");
            flint_printf("success: %d\n", success);
            flint_printf("A:\n");
            fmpz_mat_print(A);
            flint_printf("B:\n");
            fmpz_mat_print(B);
            flint_printf("X:\n");
            fmpq_mat_print(X);
            flint_printf("AX:\n");
            fmpq_mat_print(AX);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(AX_Z);
        fmpq_mat_clear(X);
        fmpq_mat_clear(AX);
    }

    /* Check singular systems */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t X;
        fmpz_mat_t A, B;
        slong n, m, bits;
        int success;

        n = 1 + n_randint(state, 40);
        m = 1 + n_randint(state, 40);
        bits = 1 + n_randint(state, 100);

        fmpz_mat_init(A, n, n);
        fmpz_mat_randrank(A, state, n_randint(state, n), bits);

        fmpz_mat_init(B, n, m);
        fmpz_mat_randtest(B, state, bits);
        fmpq_mat_init(X, n, m);

        success = fmpq_mat_solve_fmpz_mat_fraction_free(X, A, B);

        if (success != 0)
        {
            flint_printf("FAIL!\n");
            flint_printf("Expected success = 0\n");
            fmpz_mat_print(A);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpq_mat_clear(X);
    }

    TEST_FUNCTION_END(state);
}
