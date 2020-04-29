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
    
    flint_printf("add/sub....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 256);
        while (bits < UWORD(2));
        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);

        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);
        fmpz_sparse_vec_init(w);
        fmpz_sparse_vec_init(x);

        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(v, state, nnz, len, bits);

        fmpz_sparse_vec_add(w, u, v);
        fmpz_sparse_vec_sub(x, w, v);

        if (!fmpz_sparse_vec_equal(u, x, 0))
        {
            flint_printf("FAIL: u != u+v-v\n");
            flint_printf("u = "), fmpz_sparse_vec_print_pretty(u, 0, 0);
            flint_printf("v = "), fmpz_sparse_vec_print_pretty(v, 0, 0);
            flint_printf("w = "), fmpz_sparse_vec_print_pretty(w, 0, 0);
            flint_printf("x = "), fmpz_sparse_vec_print_pretty(x, 0, 0);
            abort();
        }

        fmpz_sparse_vec_add(u, u, v);
        if (!fmpz_sparse_vec_equal(u, w, 0))
        {
            flint_printf("FAIL: (u += v) != u + v\n");
            abort();
        }

        fmpz_sparse_vec_sub(u, u, v);
        if (!fmpz_sparse_vec_equal(u, x, 0))
        {
            flint_printf("FAIL: ((u += v) -= v) != u+v-v\n");
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
