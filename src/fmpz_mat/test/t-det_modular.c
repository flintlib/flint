/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_det_modular, state)
{
    fmpz_mat_t A;
    slong i, m;

    fmpz_t det1, det2;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int proved = n_randlimb(state) % 2;
        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);

        fmpz_init(det1);
        fmpz_init(det2);

        fmpz_mat_randtest(A, state, 1+n_randint(state,200));

        fmpz_mat_det_bareiss(det1, A);
        fmpz_mat_det_modular(det2, A, proved);

        if (!fmpz_equal(det1, det2))
        {
            flint_printf("FAIL:\n");
            flint_printf("different determinants!\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            flint_printf("det1: "), fmpz_print(det1), flint_printf("\n");
            flint_printf("det2: "), fmpz_print(det2), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(det1);
        fmpz_clear(det2);
        fmpz_mat_clear(A);
    }

    for (i = 0; i < 10000; i++)
    {
        int proved = n_randlimb(state) % 2;
        m = 2 + n_randint(state, 10);
        fmpz_mat_init(A, m, m);
        fmpz_init(det2);

        fmpz_mat_randrank(A, state, 1+n_randint(state, m - 1),
                                    1+n_randint(state, 10));
        fmpz_mat_randops(A, state, n_randint(state, 2*m*m + 1));

        fmpz_mat_det_modular(det2, A, proved);
        if (!fmpz_is_zero(det2))
        {
            flint_printf("FAIL:\n");
            flint_printf("expected zero determinant!\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_clear(det2);
    }

    TEST_FUNCTION_END(state);
}
