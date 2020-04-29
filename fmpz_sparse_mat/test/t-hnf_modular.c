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
#include "fmpz_sparse_mat.h"

int
main(void)
{
    slong rep, r, min_nnz, max_nnz, bits, nreps = 1000;
    fmpz_t det;
    fmpz_sparse_mat_t A, H, H2;
    FLINT_TEST_INIT(state);

    flint_printf("hnf_modular....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        bits = 1 + n_randint(state, 200);
        r = n_randint(state, 20);
        min_nnz = 0;
        max_nnz = r;

        fmpz_init(det);

        fmpz_sparse_mat_init(A, r, r);
        fmpz_sparse_mat_init(H, r, r);
        fmpz_sparse_mat_init(H2, r, r);

        do 
        {
            fmpz_sparse_mat_randtest(A, state, min_nnz, max_nnz, bits);
            fmpz_sparse_mat_det(det, A);
        } while (fmpz_is_zero(det));
        fmpz_abs(det, det);
        fmpz_mul_ui(det, det, 1 + n_randint(state, 10));

        fmpz_sparse_mat_set(H, A);
        fmpz_sparse_mat_hnf_modular(H, det);
        
        fmpz_sparse_mat_set(H2, A);
        fmpz_sparse_mat_hnf_classical(H2);
        
        if (!fmpz_sparse_mat_is_in_hnf(H))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in hnf!\n");
            fmpz_sparse_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_sparse_mat_print_pretty(H); flint_printf("\n\n");
            abort();
        }

        if (!fmpz_sparse_mat_equal(H, H2))
        {
            flint_printf("FAIL:\n");
            flint_printf("hnfs produced by different methods should be the same!\n");
            fmpz_sparse_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_sparse_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_sparse_mat_print_pretty(H2); flint_printf("\n\n");
            fmpz_print(det); flint_printf("\n\n");
            abort();
        }

        fmpz_sparse_mat_clear(H2);
        fmpz_sparse_mat_clear(H);
        fmpz_sparse_mat_clear(A);
        fmpz_clear(det);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

