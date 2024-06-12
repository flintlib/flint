/*
    Copyright (C) 2015 Alex J. Best
    Copyright (C) 2024 Joseph Grantham

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_snf_transform, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S, U, V, S2;
        slong m, n, b, d, r;
        int equal;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);
        fmpz_mat_init(U, m, m);
        fmpz_mat_init(V, n, n);
        fmpz_mat_init(S2, m, n);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, r, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        fmpz_mat_snf_transform(S, U, V, A);

        if (!fmpz_mat_is_in_snf(S))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in snf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_mul(S2, U, A);
        fmpz_mat_mul(S2, S2, V);

        if (!fmpz_mat_equal(S, S2))
        {
            flint_printf("FAIL:\n");
            flint_printf("the product of matrix with its transform matrices is not the snf!");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(V); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_snf(S2, A);

        if (!fmpz_mat_equal(S, S2))
        {
            flint_printf("FAIL:\n");
            flint_printf("the snf should agree with those produced by other methods!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(S2);
        fmpz_mat_clear(S);
        fmpz_mat_clear(U);
        fmpz_mat_clear(V);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}