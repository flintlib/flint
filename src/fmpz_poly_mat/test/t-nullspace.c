/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly_mat.h"

TEST_FUNCTION_START(fmpz_poly_mat_nullspace, state)
{
    slong i;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, N, AN;
        slong n, m, bits, deg, rank, nullity;
        float density;

        m = n_randint(state, 13);
        n = n_randint(state, 13);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(N, n, n);
        fmpz_poly_mat_init(AN, m, n);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);

        rank = fmpz_poly_mat_rank(A);
        nullity = fmpz_poly_mat_nullspace(N, A);

        if (nullity + rank != n)
        {
            flint_printf("FAIL: wrong nullity!\n");
            flint_printf("rank = %wd\n", rank);
            flint_printf("nullity = %wd\n", nullity);
            fmpz_poly_mat_print(A, "x");
            flint_printf("\n");
            fmpz_poly_mat_print(N, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        if (fmpz_poly_mat_rank(N) != nullity)
        {
            flint_printf("FAIL: wrong rank(N) != nullity!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_mul(AN, A, N);

        if (!fmpz_poly_mat_is_zero(AN))
        {
            flint_printf("FAIL: A * N != 0\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(N);
        fmpz_poly_mat_clear(AN);
    }

    TEST_FUNCTION_END(state);
}
