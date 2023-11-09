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

TEST_FUNCTION_START(fmpz_mat_snf_iliopoulos, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S, S2;
        fmpz_t mod;
        slong m, n, b, d;
        int equal;

        m = n_randint(state, 10);
        n = m;

        fmpz_init(mod);
        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);
        fmpz_mat_init(S2, m, n);

        /* sparse */
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2*m*n + 1);
        fmpz_mat_randrank(A, state, m, b);

        /* dense */
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        fmpz_mat_det(mod, A);
        fmpz_abs(mod, mod);

        fmpz_mat_snf_iliopoulos(S, A, mod);

        if (!fmpz_mat_is_in_snf(S))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in snf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_snf_iliopoulos(S2, S, mod);
        equal = fmpz_mat_equal(S, S2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("snf of a matrix in snf should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_snf_kannan_bachem(S2, S);
        equal = fmpz_mat_equal(S, S2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("snfs found by different methods should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(S2);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
        fmpz_clear(mod);
    }

    TEST_FUNCTION_END(state);
}
