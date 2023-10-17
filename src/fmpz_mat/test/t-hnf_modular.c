/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_hnf_modular, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_t det;
        fmpz_mat_t A, B, H, H2;
        slong m, n, b, c, d, i, j;
        int equal;

        n = n_randint(state, 10);
        m = n + n_randint(state, 10);

        fmpz_init(det);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        fmpz_mat_randrank(A, state, n, b);

        fmpz_mat_det(det, A);
        c = 1 + n_randint(state, 10);
        fmpz_abs(det, det);
        fmpz_mul_ui(det, det, c);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j));

        /* dense */
        d = n_randint(state, 2*m*n + 1);
        if (n_randint(state, 2))
            fmpz_mat_randops(B, state, d);

        fmpz_mat_hnf_modular(H, B, det);

        if (!fmpz_mat_is_in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_hnf_classical(H2, B);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnfs produced by different methods should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            fmpz_print(det); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_hnf_modular(H2, H, det);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("hnf of a matrix in hnf should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(H2);
        fmpz_mat_clear(H);
        fmpz_mat_clear(B);
        fmpz_mat_clear(A);
        fmpz_clear(det);
    }

    TEST_FUNCTION_END(state);
}
