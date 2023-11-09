/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_hnf_minors_transform, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, H, H2, U, UA;
        slong m, n, b, d;
        int equal;

        n = n_randint(state, 10);
        m = n + n_randint(state, 10);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);
        fmpz_mat_init(U, m, m);
        fmpz_mat_init(UA, m, n);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        fmpz_mat_randrank(A, state, n, b);

        /* dense */
        d = n_randint(state, 2*m*n + 1);
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        fmpz_mat_hnf_minors_transform(H, U, A);
        fmpz_mat_mul(UA, U, A);

        if (!fmpz_mat_is_in_hnf(H) || !fmpz_mat_equal(UA, H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_hnf_classical(H2, A);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnfs produced by different methods should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_hnf_minors(H2, H);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnf of a matrix in hnf should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(UA);
        fmpz_mat_clear(U);
        fmpz_mat_clear(H2);
        fmpz_mat_clear(H);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
