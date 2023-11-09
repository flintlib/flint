/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_solve_cramer, state)
{
    fmpz_mat_t A, X, B, AX;
    fmpz_t den;
    slong i, m, n, r;
    int success;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 4);
        n = n_randint(state, 10);

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

        success = fmpz_mat_solve_cramer(X, den, A, B);

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

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_clear(den);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 3);
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

        success = fmpz_mat_solve_cramer(X, den, A, B);

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

    TEST_FUNCTION_END(state);
}
