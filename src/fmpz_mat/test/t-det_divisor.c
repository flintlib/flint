/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_det_divisor, state)
{
    slong i;
    int result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A;
        fmpz_t det, d, q, r;
        slong m, bits;

        m = n_randint(state, 15);
        bits = 1 + n_randint(state, 50);

        fmpz_init(det);
        fmpz_init(d);
        fmpz_init(q);
        fmpz_init(r);
        fmpz_mat_init(A, m, m);

        if (i % 3 == 0 && m > 1)
        {
            /* Generate a nontrivial singular matrix */
            fmpz_mat_randrank(A, state, 1 + n_randint(state, m - 1), bits);
            fmpz_mat_randops(A, state, n_randint(state, 2*m*m + 1));
        }
        else
        {
            fmpz_mat_randtest(A, state, bits);
        }

        fmpz_mat_det_divisor(d, A);
        fmpz_mat_det_bareiss(det, A);

        if (fmpz_is_zero(det) || fmpz_is_zero(d))
        {
            result = fmpz_equal(det, d);
        }
        else
        {
            fmpz_fdiv_qr(q, r, det, d);
            result = fmpz_is_zero(r) && (fmpz_sgn(d) > 0);
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            flint_printf("det: ");  fmpz_print(det);    flint_printf("\n");
            flint_printf("d: "); fmpz_print(d); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_clear(det);
        fmpz_clear(d);
        fmpz_clear(q);
        fmpz_clear(r);
    }

    TEST_FUNCTION_END(state);
}
