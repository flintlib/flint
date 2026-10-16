/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_is_singular, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong n = n_randint(state, 12);
        slong r = n_randint(state, n + 1);
        slong bits = 1 + n_randint(state, 100);
        fmpz_mat_t A;
        int singular, reported;

        fmpz_mat_init(A, n, n);

        /* exact rank r, then random row operations */
        fmpz_mat_randrank(A, state, r, bits);
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, n_randint(state, 2 * n * n + 1));

        singular = (r < n);

        flint_fmpz_mat_force_small_primes = (int) n_randint(state, 2);
        reported = fmpz_mat_is_singular(A);
        flint_fmpz_mat_force_small_primes = 0;

        if (reported != singular)
        {
            flint_printf("FAIL\n");
            flint_printf("n = %wd, rank = %wd, reported %d\n", n, r, reported);
            fmpz_mat_print_pretty(A);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
    }


    /*
        Larger matrices, or matrices with large entries, go through the
        modular rank computation and (when singular) the exact kernel
        vector certificate rather than through a determinant.
    */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong n = 17 + n_randint(state, 8);
        slong r = n - n_randint(state, 3);
        slong bits = 1 + n_randint(state, 2) * n_randint(state, 1000);
        fmpz_mat_t A;
        int singular, reported;

        fmpz_mat_init(A, n, n);
        fmpz_mat_randrank(A, state, r, bits);
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, n_randint(state, 2 * n + 1));

        singular = (r < n);
        reported = fmpz_mat_is_singular(A);

        if (reported != singular)
        {
            flint_printf("FAIL\n");
            flint_printf("n = %wd, rank = %wd, bits = %wd, reported %d\n",
                n, r, bits, reported);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
    }

    /*
        A matrix that vanishes modulo the prime in use has rank 0 there;
        with the testing hook (p = 2) this happens for any matrix with
        even entries. With a zero column the certificate then reduces to
        checking that column; with a zero row it fails instead and the
        prime is discarded, which exercises the fallback that bounds the
        number of bad primes.
    */
    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        slong n = 17 + n_randint(state, 4);
        slong j, k;
        fmpz_mat_t A;
        int reported;

        fmpz_mat_init(A, n, n);
        fmpz_mat_randtest(A, state, 500 + n_randint(state, 500));
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                fmpz_mul_2exp(fmpz_mat_entry(A, j, k), fmpz_mat_entry(A, j, k), 1);

        if (iter % 2)
            for (j = 0; j < n; j++)
                fmpz_zero(fmpz_mat_entry(A, j, 0));   /* zero column */
        else
            for (k = 0; k < n; k++)
                fmpz_zero(fmpz_mat_entry(A, 0, k));   /* zero row */

        flint_fmpz_mat_force_small_primes = 1;
        reported = fmpz_mat_is_singular(A);
        flint_fmpz_mat_force_small_primes = 0;

        if (!reported)
        {
            flint_printf("FAIL: matrix with a zero row or column reported nonsingular\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
    }

    /*
        Very low rank: the certificate solves a small system, which is
        handled by the fraction-free algorithms rather than by lifting.
    */
    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        slong n = 17 + n_randint(state, 4);
        slong r = n_randint(state, 4);
        fmpz_mat_t A;

        fmpz_mat_init(A, n, n);
        fmpz_mat_randrank(A, state, r, 500 + n_randint(state, 500));

        if (!fmpz_mat_is_singular(A))
        {
            flint_printf("FAIL: rank %wd matrix reported nonsingular\n", r);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
