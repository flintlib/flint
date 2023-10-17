/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_det_bound, state)
{
    fmpz_mat_t A;
    slong i, m;

    fmpz_t det, bound;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);

        fmpz_init(det);
        fmpz_init(bound);

        fmpz_mat_randtest(A, state, 1+n_randint(state,200));

        fmpz_mat_det(det, A);
        fmpz_mat_det_bound(bound, A);

        if (fmpz_cmp(det, bound) > 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("bound too small!\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            flint_printf("det: "), fmpz_print(det), flint_printf("\n");
            flint_printf("bound: "), fmpz_print(bound), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(det);
        fmpz_clear(bound);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
