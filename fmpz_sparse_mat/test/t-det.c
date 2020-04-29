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
    slong rep, bits, r, nreps = 100, nzero = 0;
    fmpz_t det, d;
    fmpz_mat_t dA;
    fmpz_sparse_mat_t A;

    FLINT_TEST_INIT(state);

    flint_printf("det....");
    fflush(stdout);    

    for (rep = 0; rep < nreps; rep++)
    {
        r = n_randint(state, 75);
        bits = 1 + n_randint(state, 100);

        fmpz_sparse_mat_init(A, r, r);

        fmpz_init(det);
        fmpz_init(d);

        if (r == 0)
            fmpz_one(d);
        else 
        {
            fmpz_mat_init(dA, r, r);
            if (r > 1 && n_randint(state, 2) == 0)
            {
                fmpz_zero(d); nzero++;
                fmpz_mat_randrank(dA, state, 1+n_randint(state, r - 1), bits);
            }
            else
            {
                fmpz_randtest(d, state, 30);
                fmpz_mat_randdet(dA, state, d);        
            }
            fmpz_mat_randops(dA, state, n_randint(state, FLINT_MAX(r/4, 1)*r + 1));
            fmpz_sparse_mat_from_dense(A, dA);
            fmpz_mat_clear(dA);
        }

        fmpz_sparse_mat_det(det, A);

        if (!fmpz_equal(det, d))
        {
            flint_printf("FAIL:\n");
            flint_printf("wrong determinant!\n");
            fmpz_sparse_mat_print_pretty(A), flint_printf("\n");
            flint_printf("expected: "),  fmpz_print(d),    flint_printf("\n");
            flint_printf("ncomputed: "), fmpz_print(det), flint_printf("\n");
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
