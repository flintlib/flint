/*
    Copyright (C) 2011 Fredrik Johansson

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
    fmpz_t det, ddet;
    fmpz_mat_t dA;
    fmpz_sparse_mat_t A;

    FLINT_TEST_INIT(state);

    flint_printf("det_bareiss....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        bits = 1+n_randint(state,200);
        r = n_randint(state, 10);

        fmpz_mat_init(dA, r, r);
        fmpz_sparse_mat_init(A, r, r);
        fmpz_sparse_mat_randtest(A, state, 0, r, bits);

        fmpz_init(det);
        fmpz_init(ddet);

        fmpz_sparse_mat_det(det, A);
        fmpz_sparse_mat_to_dense(dA, A);
        fmpz_mat_det(ddet, dA);
        if (!fmpz_equal(det, ddet))
        {
            flint_printf("FAIL:\n");
            flint_printf("Incorrect determinant!\n");
            fmpz_sparse_mat_print_pretty(A), flint_printf("\n");
            flint_printf("found det: "), fmpz_print(det), flint_printf("\n");
            flint_printf("right det: "), fmpz_print(ddet), flint_printf("\n");
            abort();
        }

        fmpz_clear(det);
        fmpz_clear(ddet);
        fmpz_mat_clear(dA);
        fmpz_sparse_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
