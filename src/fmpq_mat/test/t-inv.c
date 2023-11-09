/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_inv, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C;
        fmpq_t d;

        int success1, success2;
        slong n, bits;

        n = n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, n);
        fmpq_mat_init(C, n, n);

        fmpq_init(d);

        /* XXX: replace with a randtest function */
        {
            slong k;

            for (k = 0; (k < 100) && fmpq_is_zero(d); k++)
            {
                fmpq_mat_randtest(A, state, bits);
                fmpq_mat_det(d, A);
            }
            if (fmpq_is_zero(d))
            {
                fmpq_mat_one(A);
            }
        }

        fmpq_clear(d);

        success1 = fmpq_mat_inv(B, A);
        success2 = fmpq_mat_inv(C, B);

        if (!fmpq_mat_equal(A, C) || !success1 || !success2)
        {
            flint_printf("FAIL!\n");
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

    /* Test aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B;
        fmpq_t d;

        int success1, success2;
        slong n, bits;

        n = n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, n);

        fmpq_init(d);

        /* XXX: replace with a randtest function */
        do {
            fmpq_mat_randtest(A, state, bits);
            fmpq_mat_det(d, A);
        } while (fmpq_is_zero(d));

        fmpq_clear(d);

        success1 = fmpq_mat_inv(B, A);
        success2 = fmpq_mat_inv(A, A);

        if (!fmpq_mat_equal(A, B) || !success1 || !success2)
        {
            flint_printf("FAIL!\n");
            flint_printf("A:\n");
            fmpq_mat_print(A);
            flint_printf("B:\n");
            fmpq_mat_print(B);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
    }

    /* Test singular matrices */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong n, r, b, d;
        fmpq_mat_t A, B;
        fmpz_mat_t M;
        fmpz_t den;
        int success;

        n = n_randint(state, 10);

        fmpz_init(den);

        for (r = 0; r < n; r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2*n*n + 1);

            fmpz_mat_init(M, n, n);
            fmpq_mat_init(A, n, n);
            fmpq_mat_init(B, n, n);

            fmpz_mat_randrank(M, state, r, b);

            if (i % 2 == 0)
                fmpz_mat_randops(M, state, d);

            fmpz_randtest_not_zero(den, state, b);
            fmpq_mat_set_fmpz_mat_div_fmpz(A, M, den);

            success = fmpq_mat_inv(B, A);

            if (success)
            {
                flint_printf("FAIL:\n");
                flint_printf("matrix reported as invertible:\n");
                fmpq_mat_print(A);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_clear(M);
            fmpq_mat_clear(A);
            fmpq_mat_clear(B);
        }

        fmpz_clear(den);
    }

    TEST_FUNCTION_END(state);
}
