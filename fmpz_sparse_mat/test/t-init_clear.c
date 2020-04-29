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
    slong rep, bits, r, c, i;
    fmpz_sparse_mat_t A;
    FLINT_TEST_INIT(state);
    

    flint_printf("init/clear....");
    fflush(stdout);

    for (rep = 0; rep < 100; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        r = n_randint(state, 200);
        c = n_randint(state, 200);
        fmpz_sparse_mat_init(A, r, c);

        if (!fmpz_sparse_mat_is_zero(A))
        {
            flint_printf("FAIL: A not zero!\n");
            abort();
        }
        for (i = 0; i < r; i++)
        {
            if (!fmpz_sparse_vec_is_zero(&A->rows[i]))
            {
                flint_printf("FAIL: row %wd not zero!\n", i);
                abort();
            }
        }

        fmpz_sparse_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
