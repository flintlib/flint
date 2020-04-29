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
    fmpz_sparse_mat_t A, B, C;
    FLINT_TEST_INIT(state);
    

    flint_printf("transpose....");
    fflush(stdout);

    /* Rectangular transpose, same modulus */
    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        r = n_randint(state, 20);
        c = n_randint(state, 20);
        fmpz_sparse_mat_init(A, r, c);
        fmpz_sparse_mat_init(B, c, r);
        fmpz_sparse_mat_init(C, r, c);

        fmpz_sparse_mat_randtest(A, state, 0, c, bits);
        fmpz_sparse_mat_randtest(B, state, 0, r, bits);
        fmpz_sparse_mat_randtest(C, state, 0, c, bits);

        fmpz_sparse_mat_transpose(B, A);
        fmpz_sparse_mat_transpose(C, B);
        
        if (!fmpz_sparse_mat_equal(C, A))
        {
            flint_printf("FAIL: C != A\n");
            abort();
        }

        fmpz_sparse_mat_clear(A);
        fmpz_sparse_mat_clear(B);
        fmpz_sparse_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
