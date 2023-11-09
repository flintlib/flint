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

TEST_FUNCTION_START(fmpz_mat_solve_dixon, state)
{
    fmpz_mat_t A, X, B, AX, AXm, Bm;
    fmpz_t mod;
    slong i, m, n, r;
    int success;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(Bm, m, n);
        fmpz_mat_init(X, m, n);
        fmpz_mat_init(AX, m, n);
        fmpz_mat_init(AXm, m, n);
        fmpz_init(mod);

        fmpz_mat_randrank(A, state, m, 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        success = fmpz_mat_solve_dixon(X, mod, A, B);

        fmpz_mat_set(AXm, X);

        fmpz_mat_mul(AX, A, AXm);

        fmpz_mat_scalar_mod_fmpz(AXm, AX, mod);
        fmpz_mat_scalar_mod_fmpz(Bm, B, mod);

        if (!fmpz_mat_equal(AXm, Bm) || !success)
        {
            flint_printf("FAIL:\n");
            flint_printf("AX != B!\n");
            flint_printf("A:\n"),      fmpz_mat_print_pretty(A),  flint_printf("\n");
            flint_printf("B:\n"),      fmpz_mat_print_pretty(B),  flint_printf("\n");
            flint_printf("X:\n"),      fmpz_mat_print_pretty(X),  flint_printf("\n");
            flint_printf("mod = "),    fmpz_print(mod),           flint_printf("\n");
            flint_printf("AX:\n"),     fmpz_mat_print_pretty(AX), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(Bm);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_mat_clear(AXm);
        fmpz_clear(mod);
    }

    /* Test singular systems */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 10);
        n = 1 + n_randint(state, 10);
        r = n_randint(state, m);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, m, n);
        fmpz_init(mod);

        fmpz_mat_randrank(A, state, r, 1+n_randint(state, 2)*n_randint(state, 100));
        fmpz_mat_randtest(B, state, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        if (fmpz_mat_solve_dixon(X, mod, A, B) != 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("singular system, returned nonzero\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_clear(mod);
    }

    TEST_FUNCTION_END(state);
}
