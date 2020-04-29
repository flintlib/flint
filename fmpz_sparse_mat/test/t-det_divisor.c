/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_sparse_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    slong rep, bits, r, nreps = 1000;
    fmpz_t det, d;
    fmpz_sparse_mat_t A;

    FLINT_TEST_INIT(state);

    flint_printf("det_divisor....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        r = n_randint(state, 15);
        bits = 1 + n_randint(state, 100);

        fmpz_init(det);
        fmpz_init(d);
        fmpz_sparse_mat_init(A, r, r);

        fmpz_sparse_mat_randtest(A, state, 0, r, bits);

        fmpz_sparse_mat_det_divisor(d, A);
        fmpz_sparse_mat_det_bareiss(det, A);

        if ((fmpz_is_zero(det) || fmpz_is_zero(d)) && !fmpz_equal(det, d))
        {
            flint_printf("FAIL: found divisor of matrix with det 0\n");
            fmpz_sparse_mat_print_pretty(A), flint_printf("\n");
            abort();
        }
        else if (!fmpz_divisible(det, d))
        {
            flint_printf("FAIL:\n");
            fmpz_sparse_mat_print_pretty(A), flint_printf("\n");
            flint_printf("det: ");  fmpz_print(det);    flint_printf("\n");
            flint_printf("d: "); fmpz_print(d); flint_printf("\n");
            abort();
        }

        fmpz_sparse_mat_clear(A);
        fmpz_clear(det);
        fmpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
