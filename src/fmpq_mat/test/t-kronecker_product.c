/*
    Copyright (C) 2018 Martin Raum

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_kronecker_product, state)
{
    int r, result;
    fmpq_mat_t A, B, C;
    fmpq_mat_t window1, window2;
    slong m, n, k, l, i, j;
    slong bits;

    for (r = 0; r < 100 * flint_test_multiplier(); r++)
    {
        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);
        l = n_randint(state, 10);
        if ( m && n )
        {
            i = n_randint(state, m);
            j = n_randint(state, n);
        }

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, k, l);
        fmpq_mat_init(C, m*k, n*l);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_randtest(A, state, bits);
        fmpq_mat_randtest(B, state, bits);

        fmpq_mat_kronecker_product(C, A, B);

        if ( m && n )
        {
            fmpq_mat_window_init(window1, C, 0, 0, k, l);
            fmpq_mat_window_init(window2, C, i*k, j*l, (i+1)*k, (j+1)*l);

            fmpq_mat_scalar_mul_fmpq(window1, window1, fmpq_mat_entry(A, i, j));
            fmpq_mat_scalar_mul_fmpq(window2, window2, fmpq_mat_entry(A, 0, 0));

            result = fmpq_mat_equal(window1, window2);
            if (!result)
            {
                flint_printf("FAIL:\n");
                flint_printf("A:\n");
                fmpq_mat_print(A);
                flint_printf("B:\n");
                fmpq_mat_print(B);
                flint_printf("C:\n");
                fmpq_mat_print(C);
                flint_printf("i,j: %d,%d\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpq_mat_window_clear(window1);
            fmpq_mat_window_clear(window2);
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
