/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_scalar_divexact_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t Q, A, B, Abig, Bbig;
        slong rows, cols, i, j, qbits, cbits;
        int window, alias;
        fmpz_t c;

        rows = n_randint(state, 8);
        cols = n_randint(state, 8);
        window = n_randint(state, 2);
        alias = n_randint(state, 2);

        /* divisors from one limb to many; quotients much smaller,
           comparable or larger, so that the shared inverse is computed
           to various precisions */
        cbits = n_randint(state, 4) == 0 ? 1 + n_randint(state, 3000) : 1 + n_randint(state, 300);
        qbits = n_randint(state, 3) == 0 ? 1 + n_randint(state, 3000) : 1 + n_randint(state, 200);

        fmpz_init(c);
        fmpz_randtest_not_zero(c, state, cbits);

        fmpz_mat_init(Q, rows, cols);
        fmpz_mat_randtest(Q, state, qbits);

        /* some rows of only small entries or zeros */
        for (i = 0; i < rows; i++)
        {
            if (n_randint(state, 4) == 0)
                for (j = 0; j < cols; j++)
                    fmpz_set_si(fmpz_mat_entry(Q, i, j), (slong) n_randint(state, 5) - 2);
        }

        /* A and B are either matrices or windows of larger matrices
           (with row stride larger than the number of columns) */
        if (window)
        {
            fmpz_mat_init(Abig, rows + 2, cols + 3);
            fmpz_mat_init(Bbig, rows + 1, cols + 2);
            fmpz_mat_randtest(Bbig, state, 100);
            fmpz_mat_window_init(A, Abig, 1, 2, rows + 1, cols + 2);
            fmpz_mat_window_init(B, Bbig, 1, 1, rows + 1, cols + 1);
        }
        else
        {
            fmpz_mat_init(A, rows, cols);
            fmpz_mat_init(B, rows, cols);
            fmpz_mat_randtest(B, state, 100);
        }

        fmpz_mat_scalar_mul_fmpz(A, Q, c);

        if (alias)
        {
            fmpz_mat_scalar_divexact_fmpz(A, A, c);

            if (!fmpz_mat_equal(A, Q))
                TEST_FUNCTION_FAIL("aliasing\nc = %{fmpz}\nQ = %{fmpz_mat}\nA = %{fmpz_mat}\n", c, Q, A);
        }
        else
        {
            fmpz_mat_scalar_divexact_fmpz(B, A, c);

            if (!fmpz_mat_equal(B, Q))
                TEST_FUNCTION_FAIL("c = %{fmpz}\nQ = %{fmpz_mat}\nB = %{fmpz_mat}\n", c, Q, B);
        }

        if (window)
        {
            fmpz_mat_window_clear(A);
            fmpz_mat_window_clear(B);
            fmpz_mat_clear(Abig);
            fmpz_mat_clear(Bbig);
        }
        else
        {
            fmpz_mat_clear(A);
            fmpz_mat_clear(B);
        }

        fmpz_mat_clear(Q);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
