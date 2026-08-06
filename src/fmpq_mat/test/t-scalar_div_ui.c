/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_scalar_div_ui, state)
{
    int i, result;

    /* Aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B;
        ulong c;

        slong m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);

        fmpq_mat_randtest(B, state, bits);
        c = n_randtest_not_zero(state);

        fmpq_mat_scalar_div_si(A, B, c);
        fmpq_mat_scalar_div_si(B, B, c);

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
    }

    /* (A + B) / x == A / x + B / x */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C, D;
        ulong c;

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
        c = n_randtest_not_zero(state);

        fmpq_mat_scalar_div_si(C, A, c);
        fmpq_mat_scalar_div_si(D, B, c);
        fmpq_mat_add(D, C, D);

        fmpq_mat_add(C, A, B);
        fmpq_mat_scalar_div_si(C, C, c);

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
    }

    TEST_FUNCTION_END(state);
}
