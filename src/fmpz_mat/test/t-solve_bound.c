/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_solve_bound, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, X;
        fmpz_t N, D, den;
        slong m, n, b1, b2;
        slong j, k;

        b1 = 1 + n_randint(state, 100);
        b2 = 1 + n_randint(state, 100);
        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_init(den);
        fmpz_init(N);
        fmpz_init(D);
        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, m, n);

        fmpz_mat_randrank(A, state, m, b1);
        fmpz_mat_randops(A, state, n_randint(state, m)*n_randint(state, m));
        fmpz_mat_randtest(B, state, b2);

        fmpz_mat_solve_bound(N, D, A, B);
        fmpz_mat_solve(X, den, A, B);

        if (fmpz_cmpabs(D, den) < 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("denominator bound:\n");
            fmpz_print(D);
            flint_printf("\ndenominator:\n");
            fmpz_print(den);
            flint_printf("\n");
            flint_printf("A:\n");
            fmpz_mat_print_pretty(A);
            flint_printf("B:\n");
            fmpz_mat_print_pretty(B);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        for (j = 0; j < m; j++)
        {
            for (k = 0; k < n; k++)
            {
                if (fmpz_cmpabs(N, fmpz_mat_entry(X, j, k)) < 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("numerator bound:\n");
                    fmpz_print(N);
                    flint_printf("\nnumerator:\n");
                    fmpz_print(fmpz_mat_entry(X, j, k));
                    flint_printf("\n");
                    flint_printf("A:\n");
                    fmpz_mat_print_pretty(A);
                    flint_printf("B:\n");
                    fmpz_mat_print_pretty(B);
                    flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);

        fmpz_clear(den);
        fmpz_clear(N);
        fmpz_clear(D);
    }

    TEST_FUNCTION_END(state);
}
