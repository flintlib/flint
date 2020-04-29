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
#include <limits.h>
#include "flint.h"
#include "fmpz_sparse_mat.h"

int
main(void)
{
    slong rep, bits, r, c;
    fmpz_sparse_mat_t A, B;
    fmpz_mat_t C, D;
    FLINT_TEST_INIT(state);
    

    flint_printf("conversion to/from dense matrix....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 100);
        while (bits < UWORD(2));
        r = n_randint(state, 100);
        c = n_randint(state, 100);
        fmpz_sparse_mat_init(A, r, c);
        fmpz_sparse_mat_init(B, r, c);
        fmpz_mat_init(C, r, c);
        fmpz_mat_init(D, r, c);
        
        fmpz_sparse_mat_randtest(A, state, 0, c, bits);
        fmpz_sparse_mat_to_dense(C, A);
        fmpz_sparse_mat_from_dense(B, C);
        
        if (!fmpz_sparse_mat_equal(A, B))
        {
            flint_printf("FAIL: A != B\n");
            abort();
        }

        fmpz_mat_randtest(C, state, bits);
        fmpz_sparse_mat_from_dense(A, C);
        fmpz_sparse_mat_to_dense(D, A);
        
        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: C != D\n");
            abort();
        }
        fmpz_sparse_mat_clear(A);
        fmpz_sparse_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
