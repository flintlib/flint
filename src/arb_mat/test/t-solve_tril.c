/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_solve_tril, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_mat_t A, X, B, Y;
        slong rows, cols, prec, i, j;
        int unit;

        prec = 2 + n_randint(state, 200);
        if (n_randint(state, 10) == 0)
        {
            rows = n_randint(state, 60);
            cols = n_randint(state, 60);
        }
        else
        {
            rows = n_randint(state, 10);
            cols = n_randint(state, 10);
        }
        unit = n_randint(state, 2);

        arb_mat_init(A, rows, rows);
        arb_mat_init(B, rows, cols);
        arb_mat_init(X, rows, cols);
        arb_mat_init(Y, rows, cols);

        arb_mat_randtest(A, state, prec, 10);
        arb_mat_randtest(X, state, prec, 10);
        arb_mat_randtest(Y, state, prec, 10);

        for (i = 0; i < rows; i++)
        {
            if (unit)
                arb_one(arb_mat_entry(A, i, i));
            else
                arb_set_ui(arb_mat_entry(A, i, i), 1 + n_randint(state, 100));

            for (j = i + 1; j < rows; j++)
                arb_zero(arb_mat_entry(A, i, j));
        }

        arb_mat_mul(B, A, X, prec);

        if (unit)  /* check that diagonal entries are ignored */
        {
            for (i = 0; i < rows; i++)
                arb_set_ui(arb_mat_entry(A, i, i), 1 + n_randint(state, 100));
        }

        /* Check Y = A^(-1) * (A * X) = X */
        arb_mat_solve_tril(Y, A, B, unit, prec);

        if (!arb_mat_overlaps(Y, X))
        {
            flint_printf("FAIL\n");
            flint_printf("A = \n"); arb_mat_printd(A, 10); flint_printf("\n\n");
            flint_printf("B = \n"); arb_mat_printd(B, 10); flint_printf("\n\n");
            flint_printf("X = \n"); arb_mat_printd(X, 10); flint_printf("\n\n");
            flint_printf("Y = \n"); arb_mat_printd(Y, 10); flint_printf("\n\n");
            flint_abort();
        }

        /* Check aliasing */
        arb_mat_solve_tril(B, A, B, unit, prec);
        if (!arb_mat_equal(B, Y))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("A = \n"); arb_mat_printd(A, 10); flint_printf("\n\n");
            flint_printf("B = \n"); arb_mat_printd(B, 10); flint_printf("\n\n");
            flint_printf("X = \n"); arb_mat_printd(X, 10); flint_printf("\n\n");
            flint_printf("Y = \n"); arb_mat_printd(Y, 10); flint_printf("\n\n");
            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(X);
        arb_mat_clear(Y);
    }

    TEST_FUNCTION_END(state);
}
