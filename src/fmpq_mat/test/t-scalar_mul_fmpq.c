/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2018 Martin Raum

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_scalar_mul_fmpq, state)
{
    int i, result;

    /* Aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B;
        fmpq_t x;

        slong m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_init(x);

        fmpq_mat_randtest(B, state, bits);
        fmpq_randtest(x, state, bits);

        fmpq_mat_scalar_mul_fmpq(A, B, x);
        fmpq_mat_scalar_mul_fmpq(B, B, x);

        result = fmpq_mat_equal(A, B);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n"), fmpq_mat_print(A);
            flint_printf("B:\n"), fmpq_mat_print(B);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_clear(x);
    }

    /* (A + B) * x == A * x + B * x */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C, D;
        fmpq_t x;

        slong m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_mat_init(C, m, n);
        fmpq_mat_init(D, m, n);
        fmpq_init(x);

        fmpq_mat_randtest(A, state, bits);
        fmpq_mat_randtest(B, state, bits);
        fmpq_randtest(x, state, bits);

        fmpq_mat_scalar_mul_fmpq(C, A, x);
        fmpq_mat_scalar_mul_fmpq(D, B, x);
        fmpq_mat_add(D, C, D);

        fmpq_mat_add(C, A, B);
        fmpq_mat_scalar_mul_fmpq(C, C, x);

        result = fmpq_mat_equal(C, D);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n"), fmpq_mat_print(A);
            flint_printf("B:\n"), fmpq_mat_print(B);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        fmpq_mat_clear(D);
        fmpq_clear(x);
    }

    TEST_FUNCTION_END(state);
}
