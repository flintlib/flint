/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_mul, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C, D;

        slong m, n, k, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, k);
        fmpq_mat_init(B, k, n);
        fmpq_mat_init(C, m, n);
        fmpq_mat_init(D, m, n);

        fmpq_mat_randtest(A, state, bits);
        fmpq_mat_randtest(B, state, bits);
        fmpq_mat_randtest(C, state, bits);  /* noise in output */

        fmpq_mat_mul_direct(C, A, B);
        fmpq_mat_mul_cleared(D, A, B);

        result = fmpq_mat_equal(C, D);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpq_mat_print(A);
            flint_printf("B:\n");
            fmpq_mat_print(B);
            flint_printf("C:\n");
            fmpq_mat_print(C);
            flint_printf("D:\n");
            fmpq_mat_print(D);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        fmpq_mat_clear(D);
    }

    /* Test aliasing with windows */
    {
        fmpq_mat_t A, B, A_window;

        fmpq_mat_init(A, 2, 2);
        fmpq_mat_init(B, 2, 2);

        fmpq_mat_window_init(A_window, A, 0, 0, 2, 2);

        fmpq_mat_one(A);
        fmpq_mat_one(B);
        fmpq_set_ui(fmpq_mat_entry(B, 0, 1), 1, 1);
        fmpq_set_ui(fmpq_mat_entry(B, 1, 0), 1, 1);

        fmpq_mat_mul(A_window, B, A_window);

        if (!fmpq_mat_equal(A, B))
        {
            flint_printf("FAIL: window aliasing failed\n");
            fmpq_mat_print(A); flint_printf("\n\n");
            fmpq_mat_print(B); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_window_clear(A_window);
        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
