/*
    Copyright (C) 2015 Tommy Hofmann

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
    slong rep, r, c, min_nnz, max_nnz, bits, rk, rk2, nreps = 10;
    fmpz_t mod;
    fmpz_mat_t dH;
    fmpz_sparse_mat_t A, H, H2;
    FLINT_TEST_INIT(state);

    flint_printf("howell_form....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        bits = 64 + n_randint(state, 70);
        c = n_randint(state, 10); 
        r = c + n_randint(state, 10); 
        min_nnz = 0;
        max_nnz = c;

        fmpz_mat_init(dH, r, c);
        fmpz_sparse_mat_init(A, r, c);
        fmpz_sparse_mat_init(H, r, c);
        fmpz_sparse_mat_init(H2, r, c);
        fmpz_sparse_mat_randtest(A, state, min_nnz, max_nnz, bits);

        fmpz_init(mod);
        do fmpz_randtest_not_zero(mod, state, bits);
        while (fmpz_fits_si(mod));
        fmpz_sparse_mat_set(H, A);
        fmpz_sparse_mat_to_dense(dH, A);
        rk = fmpz_sparse_mat_howell_form_mod(H, mod);
        rk2 = fmpz_mat_howell_form_mod(dH, mod);
        fmpz_sparse_mat_from_dense(H2, dH);
        if (!fmpz_sparse_mat_equal(H, H2))
        {
            flint_printf("FAIL: Howell form not equal\n");
            fmpz_sparse_mat_print_pretty(H); flint_printf("\n\n");
            fmpz_sparse_mat_print_pretty(H2); flint_printf("\n\n");
            flint_printf("Ranks: %wd, %wd\nModulus: ", rk, rk2);
            fmpz_print(mod);
            flint_printf("\n\n");
            abort();
        }

        do fmpz_randtest_unsigned(mod, state, 10); 
        while (fmpz_is_zero(mod));

        fmpz_sparse_mat_clear(A);
        fmpz_sparse_mat_clear(H);
        fmpz_sparse_mat_clear(H2);
        fmpz_mat_clear(dH);
        fmpz_clear(mod);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

