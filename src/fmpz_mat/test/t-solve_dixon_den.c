/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_solve_dixon_den, state)
{
    fmpz_mat_t A, X, B, AX;
    fmpz_t den;
    slong i, m, n, r;
    int success;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        flint_fmpz_mat_force_small_primes = n_randint(state, 2);
        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, m, n);
        fmpz_mat_init(AX, m, n);
        fmpz_init(den);

        fmpz_mat_randrank(A, state, m, 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        success = fmpz_mat_solve_dixon_den(X, den, A, B);

        fmpz_mat_mul(AX, A, X);
        fmpz_mat_scalar_divexact_fmpz(AX, AX, den);

        if (!fmpz_mat_equal(AX, B) || !success)
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


        /* the denominator must be positive, minimal (no common factor
           with the entries of X) and a divisor of det(A) */
        if (success && m > 0)
        {
            fmpz_t g, det, rem;
            slong ii, jj;
            fmpz_init(g);
            fmpz_init(det);
            fmpz_init(rem);
            fmpz_set(g, den);
            for (ii = 0; ii < m; ii++)
                for (jj = 0; jj < n; jj++)
                    fmpz_gcd(g, g, fmpz_mat_entry(X, ii, jj));
            fmpz_mat_det(det, A);
            fmpz_mod(rem, det, den);
            if (fmpz_sgn(den) <= 0 || !fmpz_is_one(g) || !fmpz_is_zero(rem))
            {
                flint_printf("FAIL:\n");
                flint_printf("denominator not positive, minimal or a divisor of det\n");
                flint_printf("A:\n"),      fmpz_mat_print_pretty(A), flint_printf("\n");
                flint_printf("X:\n"),      fmpz_mat_print_pretty(X), flint_printf("\n");
                flint_printf("den(X) = "), fmpz_print(den),          flint_printf("\n");
                flint_printf("gcd = "),    fmpz_print(g),            flint_printf("\n");
                flint_printf("det(A) = "), fmpz_print(det),          flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
            fmpz_clear(g);
            fmpz_clear(det);
            fmpz_clear(rem);
        }
        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_clear(den);
    }

    /* Test singular systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        flint_fmpz_mat_force_small_primes = n_randint(state, 2);
        m = 1 + n_randint(state, 10);
        n = 1 + n_randint(state, 10);
        r = n_randint(state, m);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, m, n);
        fmpz_mat_init(AX, m, n);
        fmpz_init(den);

        fmpz_mat_randrank(A, state, r, 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        success = fmpz_mat_solve_dixon_den(X, den, A, B);

        if (!fmpz_is_zero(den) || success)
        {
            flint_printf("FAIL:\n");
            flint_printf("singular system gave nonzero determinant\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_clear(den);
    }

    /*
        Large enough to exercise the parts of the lifting that only engage
        for long lifts: p^k digits with a Newton inverse, multimodular
        matrix-vector products through precomputed residue tables (which
        are also grown once, for cols = 8), and reconstruction with a
        precomputed inverse of the modulus. These take a few tenths of a
        second, and repeating them adds no coverage, so they run once
        regardless of the test multiplier.
    */
    {
        slong cols;

        for (cols = 1; cols <= 8; cols *= 8)
        {
            fmpz_mat_t A, B, X, AX, dB;
            fmpz_t den;
            int success;

            fmpz_mat_init(A, 24, 24);
            fmpz_mat_init(B, 24, cols);
            fmpz_mat_init(X, 24, cols);
            fmpz_mat_init(AX, 24, cols);
            fmpz_mat_init(dB, 24, cols);
            fmpz_init(den);

            fmpz_mat_randbits(A, state, 2000);
            fmpz_mat_randbits(B, state, 2000);

            success = fmpz_mat_solve_dixon_den(X, den, A, B);

            fmpz_mat_mul(AX, A, X);
            fmpz_mat_scalar_mul_fmpz(dB, B, den);

            if (!success || !fmpz_mat_equal(AX, dB))
            {
                flint_printf("FAIL:\n");
                flint_printf("large p-adic solve, cols = %wd\n", cols);
                flint_printf("success = %d, den = ", success);
                fmpz_print(den), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_clear(A);
            fmpz_mat_clear(B);
            fmpz_mat_clear(X);
            fmpz_mat_clear(AX);
            fmpz_mat_clear(dB);
            fmpz_clear(den);
        }
    }

    /*
        A solution much smaller than the matrix: the certificate is out of
        reach, so the candidate is verified by an explicit product, and for
        a small solution that product is done classically.
    */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, X, X0, AX, dB;
        fmpz_t den;
        slong m = 2 + n_randint(state, 8);
        int success;

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, 1);
        fmpz_mat_init(X, m, 1);
        fmpz_mat_init(X0, m, 1);
        fmpz_mat_init(AX, m, 1);
        fmpz_mat_init(dB, m, 1);
        fmpz_init(den);

        fmpz_mat_randbits(A, state, 500 + n_randint(state, 2000));
        fmpz_mat_randbits(X0, state, 1 + n_randint(state, 8));
        fmpz_mat_mul(B, A, X0);

        success = fmpz_mat_solve_dixon_den(X, den, A, B);

        fmpz_mat_mul(AX, A, X);
        fmpz_mat_scalar_mul_fmpz(dB, B, den);

        if (!success || !fmpz_mat_equal(AX, dB))
        {
            flint_printf("FAIL:\n");
            flint_printf("small solution of a large system\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_mat_clear(X0);
        fmpz_mat_clear(AX);
        fmpz_mat_clear(dB);
        fmpz_clear(den);
    }

    /*
        Shapes that select the remaining strategies inside the lifting: a
        larger matrix with word-size entries uses residue tables for the
        word-digit phase (rather than integer products), and a wide right
        hand side makes reconstruction attempts expensive enough that the
        probe is consulted before materialising the solution. These also
        run once.
    */
    {
        slong c;

        for (c = 0; c < 3; c++)
        {
            fmpz_mat_t A, B, X, AX, dB;
            fmpz_t den;
            slong m = (c == 2) ? 20 : 40;
            slong cols = (c == 2) ? 20 : 1;
            slong bits = (c == 0) ? 200 : ((c == 1) ? 300 : 800);
            int success;

            fmpz_mat_init(A, m, m);
            fmpz_mat_init(B, m, cols);
            fmpz_mat_init(X, m, cols);
            fmpz_mat_init(AX, m, cols);
            fmpz_mat_init(dB, m, cols);
            fmpz_init(den);

            fmpz_mat_randbits(A, state, bits);
            fmpz_mat_randbits(B, state, bits);

            success = fmpz_mat_solve_dixon_den(X, den, A, B);

            fmpz_mat_mul(AX, A, X);
            fmpz_mat_scalar_mul_fmpz(dB, B, den);

            if (!success || !fmpz_mat_equal(AX, dB))
            {
                flint_printf("FAIL:\n");
                flint_printf("solve with m = %wd, cols = %wd, bits = %wd\n", m, cols, bits);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_clear(A);
            fmpz_mat_clear(B);
            fmpz_mat_clear(X);
            fmpz_mat_clear(AX);
            fmpz_mat_clear(dB);
            fmpz_clear(den);
        }
    }

    flint_fmpz_mat_force_small_primes = 0;

    TEST_FUNCTION_END(state);
}
