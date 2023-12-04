/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_det, state)
{
    slong ix;

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
        slong m;
        fmpz_mat_t A;
        fmpz_t det, result;
        ulong is_small, is_singular;

        is_singular = n_randint(state, 2);
        is_small = n_randint(state, 20);

        m = n_randint(state, 10);
        if (is_singular)
            m += 2;
        if (!is_small)
            m += 60; /* Generate last conditional in fmpz_det source code */

        fmpz_mat_init(A, m, m);
        fmpz_init(det);
        fmpz_init(result);

        if (is_singular)
        {
            fmpz_zero(det);

            fmpz_mat_randrank(A, state, 1 + n_randint(state, m - 1), 1 + n_randint(state, 10));
            fmpz_mat_randops(A, state, n_randint(state, 2 * m * m + 1));
        }
        else
        {
            if (m)
                fmpz_randtest(det, state, 30);
            else
                fmpz_set_ui(det, UWORD(1));

            fmpz_mat_randdet(A, state, det);
            fmpz_mat_randops(A, state, n_randint(state, 2*m*m + 1));
        }

        fmpz_mat_det(result, A);

        if (!fmpz_equal(det, result))
        {
            flint_printf("FAIL:\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            flint_printf("expected: "),  fmpz_print(det),    flint_printf("\n");
            flint_printf("ncomputed: "), fmpz_print(result), flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_clear(det);
        fmpz_clear(result);
    }

    TEST_FUNCTION_END(state);
}
