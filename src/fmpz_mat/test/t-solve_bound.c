/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_solve_bound, state)
{
    slong i;

    {
        /* Regression test: the numerator bound must be built from the
           column Hadamard product of A. Here the row Hadamard product is
           sqrt(70) = 8.37 and the true numerators reach 90, so using
           min(row, col) for N (which is correct for D) would make the
           modular solvers, which never reconstruct beyond N, fail to
           terminate correctly. Random matrices with balanced row norms do
           not exhibit this. */
        static const slong Aentries[9] = { -1, -3, -5, 0, 1, 0, 0, 1, -1 };
        static const slong Bentries[3] = { -1, -8, 5 };
        static const slong Xentries[3] = { 90, -8, -13 };
        fmpz_mat_t A, B, X;
        fmpz_t N, D, den;
        slong j, k;

        fmpz_mat_init(A, 3, 3);
        fmpz_mat_init(B, 3, 1);
        fmpz_mat_init(X, 3, 1);
        fmpz_init(N);
        fmpz_init(D);
        fmpz_init(den);

        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
                fmpz_set_si(fmpz_mat_entry(A, j, k), Aentries[3 * j + k]);
            fmpz_set_si(fmpz_mat_entry(B, j, 0), Bentries[j]);
        }

        fmpz_mat_solve_bound(N, D, A, B);

        for (j = 0; j < 3; j++)
        {
            if (fmpz_cmp_si(N, FLINT_ABS(Xentries[j])) < 0)
            {
                flint_printf("FAIL: numerator bound too small\n");
                flint_printf("N = "), fmpz_print(N), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* the solvers must also agree with the known solution */
        for (k = 0; k < 3; k++)
        {
            int success;

            if (k == 0)
                success = fmpz_mat_solve_dixon_den(X, den, A, B);
            else if (k == 1)
                success = fmpz_mat_solve_multi_mod_den(X, den, A, B);
            else
                success = fmpz_mat_solve(X, den, A, B);

            if (!success || !fmpz_is_one(den))
            {
                flint_printf("FAIL: solve failed on the regression input\n");
                fflush(stdout);
                flint_abort();
            }

            for (j = 0; j < 3; j++)
            {
                if (!fmpz_equal_si(fmpz_mat_entry(X, j, 0), Xentries[j]))
                {
                    flint_printf("FAIL: wrong solution on the regression input\n");
                    fmpz_mat_print_pretty(X);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_clear(N);
        fmpz_clear(D);
        fmpz_clear(den);
    }

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

        if (i % 2)
        {
            fmpz_mat_randrank(A, state, m, b1);
            fmpz_mat_randops(A, state, n_randint(state, m)*n_randint(state, m));
            fmpz_mat_randtest(B, state, b2);
        }
        else
        {
            /* One dense row and sparse rows of unit entries elsewhere: the
               row Hadamard product is then much smaller than the column
               one, and the numerators can exceed the row bound (they are
               only bounded by sqrt(m) times it). Generators that keep the
               row norms balanced never exercise this, which is why the
               numerator bound must be built from the column Hadamard
               product alone. */
            m = FLINT_MIN(m, 5);
            b2 = 1 + n_randint(state, 5);
            fmpz_mat_clear(A);
            fmpz_mat_clear(B);
            fmpz_mat_clear(X);
            fmpz_mat_init(A, m, m);
            fmpz_mat_init(B, m, n);
            fmpz_mat_init(X, m, n);

            for (j = 1; j < m; j++)
            {
                fmpz_set_si(fmpz_mat_entry(A, j, j), 1 - 2 * (slong) n_randint(state, 2));
                if (n_randint(state, 2))
                    fmpz_set_si(fmpz_mat_entry(A, j, n_randint(state, m)),
                                1 - 2 * (slong) n_randint(state, 2));
            }
            for (j = 0; j < m; j++)
                fmpz_randtest(fmpz_mat_entry(A, 0, j), state, 1 + n_randint(state, 4));

            fmpz_mat_randtest(B, state, b2);
        }

        fmpz_mat_solve_bound(N, D, A, B);

        /* the sparse generator below can produce a singular matrix (a zero
           bound means a zero row or column), for which the bounds say
           nothing */
        if (fmpz_is_zero(D) || !fmpz_mat_solve(X, den, A, B))
            goto next_iter;

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

next_iter:
        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);

        fmpz_clear(den);
        fmpz_clear(N);
        fmpz_clear(D);
    }

    TEST_FUNCTION_END(state);
}
