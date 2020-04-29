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
#include <limits.h>
#include "flint.h"
#include "fmpz_sparse_mat.h"

int
main(void)
{
    slong rep, bits, r, c, nreps = 1000;
    fmpz_t a, cc;
    fmpz_sparse_mat_t A, B, C, D;
    FLINT_TEST_INIT(state);
    
    flint_printf("scalar_mul....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        r = n_randint(state, 200);
        c = n_randint(state, 200);
        fmpz_init(a);
        fmpz_init(cc);
        fmpz_randtest(cc, state, bits);
        fmpz_sparse_mat_init(A, r, c);
        fmpz_sparse_mat_init(B, r, c);
        fmpz_sparse_mat_init(C, r, c);
        fmpz_sparse_mat_init(D, r, c);

        fmpz_sparse_mat_randtest(A, state, 0, c, bits);

        fmpz_sparse_mat_scalar_mul_fmpz(B, A, a);
        fmpz_one(cc);
        fmpz_sub(cc, a, cc);
        fmpz_sparse_mat_scalar_mul_fmpz(C, A, cc);

        /* c*A - (c-1)*A == A */
        fmpz_sparse_mat_sub(D, B, C);

        if (!fmpz_sparse_mat_equal(A, D))
        {
            flint_printf("FAIL\n");
            abort();
        }

         /* Aliasing */
        fmpz_sparse_mat_scalar_mul_fmpz(C, A, a);
        fmpz_sparse_mat_scalar_mul_fmpz(A, A, a);

        if (!fmpz_sparse_mat_equal(A, C))
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_sparse_mat_clear(A);
        fmpz_sparse_mat_clear(B);
        fmpz_sparse_mat_clear(C);
        fmpz_sparse_mat_clear(D);
        fmpz_clear(a);
        fmpz_clear(cc);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
