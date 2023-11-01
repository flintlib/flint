/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "d_mat.h"

TEST_FUNCTION_START(d_mat_transpose, state)
{
    slong m, n, rep;

    /* Rectangular transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        d_mat_t A, B, C;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        d_mat_init(A, m, n);
        d_mat_init(B, n, m);
        d_mat_init(C, m, n);

        d_mat_randtest(A, state, n_randint(state, 100), n_randint(state, 100));
        d_mat_randtest(B, state, n_randint(state, 100), n_randint(state, 100));

        d_mat_transpose(B, A);
        d_mat_transpose(C, B);

        if (!d_mat_equal(C, A))
        {
            flint_printf("FAIL: C != A\n");
            flint_printf("C:\n");
            d_mat_print(C);
            flint_printf("A:\n");
            d_mat_print(A);
            fflush(stdout);
            flint_abort();
        }

        d_mat_clear(A);
        d_mat_clear(B);
        d_mat_clear(C);
    }

    /* Self-transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        d_mat_t A, B;

        m = n_randint(state, 20);

        d_mat_init(A, m, m);
        d_mat_init(B, m, m);

        d_mat_randtest(A, state, n_randint(state, 100), n_randint(state, 100));
        d_mat_set(B, A);
        d_mat_transpose(B, B);
        d_mat_transpose(B, B);

        if (!d_mat_equal(B, A))
        {
            flint_printf("FAIL: B != A\n");
            flint_printf("B:\n");
            d_mat_print(B);
            flint_printf("A:\n");
            d_mat_print(A);
            fflush(stdout);
            flint_abort();
        }

        d_mat_clear(A);
        d_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
