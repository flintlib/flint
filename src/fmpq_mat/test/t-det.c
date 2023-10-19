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

TEST_FUNCTION_START(fmpq_mat_det, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C;
        fmpq_t a, b, ab, c;

        slong n, bits;

        n = n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, n);
        fmpq_mat_init(C, n, n);

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(ab);
        fmpq_init(c);

        fmpq_mat_randtest(A, state, bits);
        fmpq_mat_randtest(B, state, bits);
        fmpq_mat_mul(C, A, B);

        fmpq_mat_det(a, A);
        fmpq_mat_det(b, B);
        fmpq_mat_det(c, C);

        fmpq_mul(ab, a, b);

        if (!fmpq_equal(ab, c))
        {
            flint_printf("FAIL!\n");
            flint_printf("A:\n");
            fmpq_mat_print(A);
            flint_printf("B:\n");
            fmpq_mat_print(B);
            flint_printf("C:\n");
            fmpq_mat_print(C);
            flint_printf("\ndet(A):\n");
            fmpq_print(a);
            flint_printf("\ndet(B):\n");
            fmpq_print(b);
            flint_printf("\ndet(C):\n");
            fmpq_print(c);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(ab);
        fmpq_clear(c);

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
