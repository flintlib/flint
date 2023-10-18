/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_can_solve_multi_mod_den, state)
{
    fmpz_mat_t A, X, B, AX;
    fmpz_t den;
    slong i, m, n, k;
    int success;

    /* test random systems (likely not soluble) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        k = n_randint(state, 20);

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, k, n);
        fmpz_mat_init(AX, m, n);
        fmpz_init(den);

        fmpz_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1), 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        success = fmpz_mat_can_solve_multi_mod_den(X, den, A, B);

        if (success)
        {
            fmpz_mat_mul(AX, A, X);
            fmpz_mat_scalar_divexact_fmpz(AX, AX, den);
        }

        if (success && !fmpz_mat_equal(AX, B))
        {
            flint_printf("FAIL:\n");
            flint_printf("AX != B!\n");
            flint_printf("A:\n"),      fmpz_mat_print_pretty(A),  flint_printf("\n");
            flint_printf("B:\n"),      fmpz_mat_print_pretty(B),  flint_printf("\n");
            flint_printf("X:\n"),      fmpz_mat_print_pretty(X),  flint_printf("\n");
            flint_printf("den(X) = "), fmpz_print(den),           flint_printf("\n");
            flint_printf("AX:\n"),     fmpz_mat_print_pretty(AX), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_clear(den);
    }

    /* test random soluble systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        k = n_randint(state, 20);

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, k, n);
        fmpz_mat_init(AX, m, n);
        fmpz_init(den);

        fmpz_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1), 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(X, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        fmpz_mat_mul(B, A, X);

        fmpz_randtest_not_zero(den, state, 20);
        fmpz_mat_scalar_mul_fmpz(A, A, den);

        success = fmpz_mat_can_solve_multi_mod_den(X, den, A, B);

        if (success)
        {
            fmpz_mat_mul(AX, A, X);
            fmpz_mat_scalar_divexact_fmpz(AX, AX, den);
        }

        if (!success || !fmpz_mat_equal(AX, B))
        {
            flint_printf("FAIL:\n");
            flint_printf("AX != B!\n");
            flint_printf("A:\n"),      fmpz_mat_print_pretty(A),  flint_printf("\n");
            flint_printf("B:\n"),      fmpz_mat_print_pretty(B),  flint_printf("\n");
            flint_printf("X:\n"),      fmpz_mat_print_pretty(X),  flint_printf("\n");
            flint_printf("den(X) = "), fmpz_print(den),           flint_printf("\n");
            flint_printf("AX:\n"),     fmpz_mat_print_pretty(AX), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_clear(den);
    }

    TEST_FUNCTION_END(state);
}
