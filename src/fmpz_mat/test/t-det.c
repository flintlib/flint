/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    fmpz_mat_t A;
    slong i, m;

    fmpz_t det, result;

    FLINT_TEST_INIT(state);

    flint_printf("det....");
    fflush(stdout);    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);

        fmpz_init(det);
        fmpz_init(result);

        if (m)
            fmpz_randtest(det, state, 30);
        else
            fmpz_set_ui(det, UWORD(1));

        fmpz_mat_randdet(A, state, det);
        fmpz_mat_randops(A, state, n_randint(state, 2*m*m + 1));

        fmpz_mat_det(result, A);

        if (!fmpz_equal(det, result))
        {
            flint_printf("FAIL:\n");
            flint_printf("wrong determinant!\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            flint_printf("expected: "),  fmpz_print(det),    flint_printf("\n");
            flint_printf("ncomputed: "), fmpz_print(result), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_clear(det);
        fmpz_clear(result);
    }

    /* Generate nontrivial singular matrices */
    for (i = 0; i < 10000; i++)
    {
        m = 2 + n_randint(state, 10);
        fmpz_mat_init(A, m, m);
        fmpz_init(det);

        fmpz_mat_randrank(A, state, 1+n_randint(state, m - 1), 1+n_randint(state, 10));
        fmpz_mat_randops(A, state, n_randint(state, 2*m*m + 1));

        fmpz_mat_det(det, A);
        if (*det)
        {
            flint_printf("FAIL:\n");
            flint_printf("expected zero determinant!\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_clear(det);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
