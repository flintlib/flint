/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_pascal, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_mat_t A;
        fmpz_t t;
        slong n, m, i, j, prec;

        n = n_randint(state, 10);
        m = n_randint(state, 10);
        prec = 2 + n_randint(state, 200);

        fmpz_init(t);
        arb_mat_init(A, n, m);
        arb_mat_randtest(A, state, 100, 10);

        arb_mat_pascal(A, 0, prec);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                fmpz_bin_uiui(t, i + j, i);
                if (!arb_contains_fmpz(arb_mat_entry(A, i, j), t))
                {
                    flint_printf("FAIL: containment (0)\n");
                    flint_abort();
                }
            }
        }

        arb_mat_pascal(A, 1, prec);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                fmpz_bin_uiui(t, j, i);
                if (!arb_contains_fmpz(arb_mat_entry(A, i, j), t))
                {
                    flint_printf("FAIL: containment (1)\n");
                    flint_abort();
                }
            }
        }

        arb_mat_pascal(A, -1, prec);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                fmpz_bin_uiui(t, i, j);
                if (!arb_contains_fmpz(arb_mat_entry(A, i, j), t))
                {
                    flint_printf("FAIL: containment (-1)\n");
                    flint_abort();
                }
            }
        }

        arb_mat_clear(A);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
