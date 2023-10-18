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

TEST_FUNCTION_START(fmpz_mat_hnf_transform, state)
{
    slong iter;

    /* Random rectangular matrices */
    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, H, H2, U;
        fmpz_t det;
        slong m, n, b, d, r;
        int equal;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);
        fmpz_mat_init(U, m, m);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, r, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        fmpz_mat_hnf_transform(H, U, A);

        if (!fmpz_mat_is_in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_init(det);

        fmpz_mat_det(det, U);
        if (!fmpz_is_pm1(det))
        {
            flint_printf("FAIL:\n");
            flint_printf("transformation matrices should have determinant +-1, U does not!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(det);

        fmpz_mat_mul(H2, U, A);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("multiplying by the transformation matrix should give the same HNF!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_hnf(H2, A);
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

        fmpz_mat_clear(U);
        fmpz_mat_clear(H2);
        fmpz_mat_clear(H);
        fmpz_mat_clear(A);
    }

    /* Now matrices with full column rank */
    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, H, H2, U;
        fmpz_t det;
        slong m, n, b, d, r;
        int equal;

        n = n_randint(state, 10);
        m = n + n_randint(state, 10);
        r = n;

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(H, m, n);
        fmpz_mat_init(H2, m, n);
        fmpz_mat_init(U, m, m);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, r, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        fmpz_mat_hnf_transform(H, U, A);

        if (!fmpz_mat_is_in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_init(det);

        fmpz_mat_det(det, U);
        if (!fmpz_is_pm1(det))
        {
            flint_printf("FAIL:\n");
            flint_printf("transformation matrices should have determinant +-1, U does not!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(det);

        fmpz_mat_mul(H2, U, A);
        equal = fmpz_mat_equal(H, H2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("multiplying by the transformation matrix should give the same HNF!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_mat_print_pretty(H2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_hnf(H2, A);
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

        fmpz_mat_clear(U);
        fmpz_mat_clear(H2);
        fmpz_mat_clear(H);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
