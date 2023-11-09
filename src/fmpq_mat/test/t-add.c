/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_add, state)
{
    int i, result;

    /* Aliasing, B = B + C */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C;

        slong m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_mat_init(C, m, n);

        fmpq_mat_randtest(B, state, bits);
        fmpq_mat_randtest(C, state, bits);

        fmpq_mat_add(A, B, C);
        fmpq_mat_add(B, B, C);

        result = fmpq_mat_equal(A, B);
        if (!result)
        {
            flint_printf("FAIL (B = B + C):\n");
            flint_printf("A:\n");
            fmpq_mat_print(A);
            flint_printf("B:\n");
            fmpq_mat_print(B);
            flint_printf("C:\n");
            fmpq_mat_print(C);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
    }

    /* Aliasing, C = B + C */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C;

        slong m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_mat_init(C, m, n);

        fmpq_mat_randtest(B, state, bits);
        fmpq_mat_randtest(C, state, bits);

        fmpq_mat_add(A, B, C);
        fmpq_mat_add(C, B, C);

        result = fmpq_mat_equal(A, C);
        if (!result)
        {
            flint_printf("FAIL (C = B + C):\n");
            flint_printf("A:\n");
            fmpq_mat_print(A);
            flint_printf("B:\n");
            fmpq_mat_print(B);
            flint_printf("C:\n");
            fmpq_mat_print(C);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
    }

    /* A + B == B + A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C, D;

        slong m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_mat_init(C, m, n);
        fmpq_mat_init(D, m, n);

        fmpq_mat_randtest(A, state, bits);
        fmpq_mat_randtest(B, state, bits);

        fmpq_mat_add(C, A, B);
        fmpq_mat_add(D, B, A);

        result = fmpq_mat_equal(C, D);
        if (!result)
        {
            flint_printf("FAIL (A + B == B + A):\n");
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

    TEST_FUNCTION_END(state);
}
