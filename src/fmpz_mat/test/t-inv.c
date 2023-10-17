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
#include "fmpz_vec.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_inv, state)
{
    fmpz_mat_t A, B, C, I;
    fmpz_t den;
    slong i, j, m, r;

    {
        fmpz_t d;
        fmpz_mat_t A, B, C;

        fmpz_mat_init(A, 1, 1);
        fmpz_one(fmpz_mat_entry(A, 0, 0));

        fmpz_mat_window_init(B, A, 0, 0, 1, 1);

        fmpz_mat_init(C, 1, 1);
        fmpz_init(d);

        fmpz_mat_inv(C, d, B);

        fmpz_clear(d);
        fmpz_mat_clear(C);
        fmpz_mat_window_clear(B);
        fmpz_mat_clear(A);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);
        fmpz_mat_init(C, m, m);
        fmpz_mat_init(I, m, m);
        fmpz_init(den);

        for (j = 0; j < m; j++)
            fmpz_set_ui(&I->rows[j][j], UWORD(1));

        /* Verify that A * A^-1 = I for random matrices */

        fmpz_mat_randrank(A, state, m, 1+n_randint(state, 2)*n_randint(state, 100));
        /* Dense or sparse? */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        fmpz_mat_inv(B, den, A);
        fmpz_mat_mul(C, A, B);

        _fmpz_vec_scalar_divexact_fmpz(C->entries, C->entries, m*m, den);

        if (!fmpz_mat_equal(C, I))
        {
            flint_printf("FAIL:\n");
            flint_printf("A * A^-1 != I!\n");
            flint_printf("A:\n"),         fmpz_mat_print_pretty(A), flint_printf("\n");
            flint_printf("A^-1:\n"),      fmpz_mat_print_pretty(B), flint_printf("\n");
            flint_printf("den(A^-1) = "), fmpz_print(den), flint_printf("\n");
            flint_printf("A * A^-1:\n"),  fmpz_mat_print_pretty(C), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* Test aliasing */
        fmpz_mat_set(C, A);
        fmpz_mat_inv(A, den, A);
        fmpz_mat_mul(B, A, C);
        _fmpz_vec_scalar_divexact_fmpz(B->entries, B->entries, m*m, den);

        if (!fmpz_mat_equal(B, I))
        {
            flint_printf("FAIL:\n");
            flint_printf("aliasing failed!\n");
            fmpz_mat_print(C); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(I);
        fmpz_clear(den);
    }

    /* Test singular matrices */
    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 10);
        r = n_randint(state, m);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);
        fmpz_init(den);

        fmpz_mat_randrank(A, state, r, 1+n_randint(state, 2)*n_randint(state, 100));

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, 1+n_randint(state, 1 + m*m));

        fmpz_mat_inv(B, den, A);
        if (!fmpz_is_zero(den))
        {
            flint_printf("FAIL:\n");
            flint_printf("singular system gave nonzero denominator\n");
            fflush(stdout);
            flint_abort();
        }

        /* Aliasing */
        fmpz_mat_inv(A, den, A);
        if (!fmpz_is_zero(den))
        {
            flint_printf("FAIL:\n");
            flint_printf("singular system gave nonzero denominator\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_clear(den);
    }

    TEST_FUNCTION_END(state);
}
