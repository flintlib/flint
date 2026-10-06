/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_hilbert, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_mat_t A;
        arb_t t, u;
        fmpq_t q;
        slong n, m, i, j, prec;

        n = n_randint(state, 12);
        m = n_randint(state, 12);
        prec = 2 + n_randint(state, 200);

        arb_init(t);
        arb_init(u);
        fmpq_init(q);
        arb_mat_init(A, n, m);
        arb_mat_randtest(A, state, 100, 10);

        arb_mat_hilbert(A, prec);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                fmpq_set_si(q, 1, i + j + 1);

                if (!arb_contains_fmpq(arb_mat_entry(A, i, j), q))
                    TEST_FUNCTION_FAIL("containment: n = %wd, m = %wd, i = %wd, j = %wd, prec = %wd\n%{arb}\n",
                        n, m, i, j, prec, arb_mat_entry(A, i, j));

                /* reference: what the old implementation computed */
                arb_one(t);
                arb_div_ui(t, t, i + j + 1, prec);

                arb_set_ui(u, i + j + 1);
                arb_ui_div(u, 1, u, prec);

                if (!arb_equal(arb_mat_entry(A, i, j), t) || !arb_overlaps(arb_mat_entry(A, i, j), u))
                    TEST_FUNCTION_FAIL("equality: n = %wd, m = %wd, i = %wd, j = %wd, prec = %wd\n%{arb}\n%{arb}\n%{arb}\n",
                        n, m, i, j, prec, arb_mat_entry(A, i, j), t, u);
            }
        }

        arb_mat_clear(A);
        arb_clear(t);
        arb_clear(u);
        fmpq_clear(q);
    }

    TEST_FUNCTION_END(state);
}
