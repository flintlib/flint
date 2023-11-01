/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_nullspace, state)
{
    fmpz_mat_t A, B, ker;
    slong i, m, n, b, d, r, nullity, nulrank;

    /* small dimension */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2*m*n + 1);

            fmpz_mat_init(A, m, n);
            fmpz_mat_init(ker, n, n);
            fmpz_mat_init(B, m, n);

            fmpz_mat_randrank(A, state, r, b);
            /* Densify */
            if (n_randlimb(state) % 2)
                fmpz_mat_randops(A, state, d);

            nullity = fmpz_mat_nullspace(ker, A);
            nulrank = fmpz_mat_rank(ker);

            if (nullity != nulrank)
            {
                flint_printf("FAIL:\n");
                flint_printf("rank(ker) != nullity!\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (nullity + r != n)
            {
                flint_printf("FAIL:\n");
                flint_printf("nullity + rank != n\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_mul(B, A, ker);

            if (fmpz_mat_rank(B) != 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("A * ker != 0\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_clear(A);
            fmpz_mat_clear(ker);
            fmpz_mat_clear(B);
        }
    }

    /* larger dimension */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        m = 25 + n_randint(state, 10);
        n = 25 + n_randint(state, 10);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2*m*n + 1);

            fmpz_mat_init(A, m, n);
            fmpz_mat_init(ker, n, n);
            fmpz_mat_init(B, m, n);

            fmpz_mat_randrank(A, state, r, b);
            /* Densify */
            if (n_randlimb(state) % 2)
                fmpz_mat_randops(A, state, d);

            nullity = fmpz_mat_nullspace(ker, A);
            nulrank = fmpz_mat_rank(ker);

            if (nullity != nulrank)
            {
                flint_printf("FAIL:\n");
                flint_printf("rank(ker) != nullity!\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (nullity + r != n)
            {
                flint_printf("FAIL:\n");
                flint_printf("nullity + rank != n\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_mul(B, A, ker);

            if (fmpz_mat_rank(B) != 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("A * ker != 0\n");
                fmpz_mat_print_pretty(A);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_clear(A);
            fmpz_mat_clear(ker);
            fmpz_mat_clear(B);
        }
    }

    TEST_FUNCTION_END(state);
}
