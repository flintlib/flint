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
    fmpz_t det, d;
    fmpz_sparse_mat_t A;

    FLINT_TEST_INIT(state);

    flint_printf("det_modular....");
    fflush(stdout);    

    for (rep = 0; rep < nreps; rep++)
    {
        int proved = n_randlimb(state) % 2;
        bits = 1+n_randint(state,200);
        r = n_randint(state, 20);

        fmpz_sparse_mat_init(A, r, r);

        fmpz_init(det);
        fmpz_init(d);

        fmpz_sparse_mat_randtest(A, state, 0, r, bits);

        fmpz_sparse_mat_det_bareiss(d, A);
        fmpz_sparse_mat_det_modular(det, A, proved);

        if (!fmpz_equal(d, det))
        {
            flint_printf("FAIL:\n");
            flint_printf("different determinants!\n");
            fmpz_sparse_mat_print_pretty(A), flint_printf("\n");
            flint_printf("det_bareiss: "), fmpz_print(d), flint_printf("\n");
            flint_printf("det_modular: "), fmpz_print(det), flint_printf("\n");
            abort();
        }

        fmpz_clear(d);
        fmpz_clear(det);
        fmpz_sparse_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
