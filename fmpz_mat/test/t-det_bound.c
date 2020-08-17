/*
    Copyright (C) 2011 Fredrik Johansson

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

    fmpz_t det, bound;

    FLINT_TEST_INIT(state);

    flint_printf("det_bound....");
    fflush(stdout);

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
            abort();
        }

        fmpz_clear(det);
        fmpz_clear(bound);
        fmpz_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
