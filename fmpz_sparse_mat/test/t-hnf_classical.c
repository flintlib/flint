/*
    Copyright (C) 2014 Alex J. Best

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
#include "fmpz_sparse_mat.h"

int
main(void)
{
    slong rep, r, c, min_nnz, max_nnz, bits, rk, nreps = 1000;
    fmpz_mat_t dA, dH;
    fmpz_sparse_mat_t A, H, H2;
    FLINT_TEST_INIT(state);

    flint_printf("hnf_classical....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        bits = 2 + n_randint(state, 100);
        r = n_randint(state, 10);
        c = n_randint(state, 10);
        min_nnz = 0;
        max_nnz = c;

        fmpz_mat_init(dA, r, c);
        fmpz_mat_init(dH, r, c);
        fmpz_sparse_mat_init(A, r, c);
        fmpz_sparse_mat_init(H, r, c);
        fmpz_sparse_mat_init(H2, r, c);
        fmpz_sparse_mat_randtest(A, state, min_nnz, max_nnz, bits);
        fmpz_sparse_mat_set(H, A);
        rk = fmpz_sparse_mat_hnf_classical(H);
        fmpz_sparse_mat_to_dense(dA, A);
        fmpz_mat_hnf_classical(dH, dA);
        fmpz_sparse_mat_from_dense(H2, dH);

        if (!fmpz_sparse_mat_equal(H, H2))
        {
            flint_printf("FAIL: matrix not in hnf!\n");
            fmpz_sparse_mat_print_pretty(A); flint_printf("\n\n");
            flint_printf("Obtained: \n");
            fmpz_sparse_mat_print_pretty(H); flint_printf("\n\n");
            flint_printf("Should be: \n");
            fmpz_sparse_mat_print_pretty(H2); flint_printf("\n\n");
            flint_printf("Found rank %wd\n", rk);
            abort();
        }

        fmpz_sparse_mat_clear(H2);
        fmpz_sparse_mat_clear(H);
        fmpz_sparse_mat_clear(A);
        fmpz_mat_clear(dA);
        fmpz_mat_clear(dH);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

