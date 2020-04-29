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
#include <gmp.h>
#include "flint.h"
#include "fmpz_sparse_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    slong rep, bits, len, nnz;
    fmpz_sparse_vec_t u, v, w, x;
    FLINT_TEST_INIT(state);
    

    flint_printf("neg....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);

        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);
        fmpz_sparse_vec_init(w);
        fmpz_sparse_vec_init(x);

        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(v, state, nnz, len, bits);

        fmpz_sparse_vec_sub(w, u, v);
        fmpz_sparse_vec_neg(v, v);
        fmpz_sparse_vec_add(x, u, v);

        if (!fmpz_sparse_vec_equal(w, x, 0))
        {
            flint_printf("FAIL: u - v != u + (-v)\n");
            abort();
        }

        fmpz_sparse_vec_neg(w, u);
        fmpz_sparse_vec_neg(u, u);

        if (!fmpz_sparse_vec_equal(w, w, 0))
        {
            flint_printf("FAIL\n");
            abort();
        }
        fmpz_sparse_vec_clear(u);
        fmpz_sparse_vec_clear(v);
        fmpz_sparse_vec_clear(w);
        fmpz_sparse_vec_clear(x);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
