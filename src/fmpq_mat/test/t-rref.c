/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_rref, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong m, n, r, rank, b, d;
        fmpq_mat_t A, B, C;
        fmpz_mat_t M;
        fmpz_t den;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_init(den);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2*m*n + 1);

            fmpz_mat_init(M, m, n);
            fmpq_mat_init(A, m, n);
            fmpq_mat_init(B, m, n);
            fmpq_mat_init(C, m, n);

            fmpz_mat_randrank(M, state, r, b);

            if (i % 2 == 0)
                fmpz_mat_randops(M, state, d);

            fmpz_randtest_not_zero(den, state, b);
            fmpq_mat_set_fmpz_mat_div_fmpz(A, M, den);

            rank = fmpq_mat_rref_classical(B, A);
            if (r != rank)
            {
                flint_printf("FAIL:\n");
                flint_printf("fmpq_mat_rref_classical: wrong rank!\n");
                fmpq_mat_print(A);
                flint_printf("\nrank: %wd computed: %wd\n", r, rank);
                fflush(stdout);
                flint_abort();
            }

            rank = fmpq_mat_rref_fraction_free(C, A);
            if (r != rank)
            {
                flint_printf("FAIL:\n");
                flint_printf("fmpq_mat_rref_fraction_free: wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }

            if (!fmpq_mat_equal(B, C))
            {
                flint_printf("FAIL:\n");
                flint_printf("different results!\n");
                flint_printf("A:\n");
                fmpq_mat_print(A);
                flint_printf("\nB:\n");
                fmpq_mat_print(B);
                flint_printf("\nC:\n");
                fmpq_mat_print(C);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_clear(M);
            fmpq_mat_clear(A);
            fmpq_mat_clear(B);
            fmpq_mat_clear(C);
        }

        fmpz_clear(den);
    }

    TEST_FUNCTION_END(state);
}
