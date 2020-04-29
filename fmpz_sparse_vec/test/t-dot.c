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
    fmpz_t a, b;
    fmpz_sparse_vec_t u, v;
    fmpz *w, *x;
    FLINT_TEST_INIT(state);
    

    flint_printf("dot product....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);

        fmpz_init(a);
        fmpz_init(b);
        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);
        w = _fmpz_vec_init(len);
        x = _fmpz_vec_init(len);

        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(v, state, nnz, len, bits);
        fmpz_sparse_vec_to_dense(w, u, len);
        fmpz_sparse_vec_to_dense(x, v, len);

        fmpz_sparse_vec_dot(a, u, v);
        _fmpz_vec_dot(b, w, x, len);

        if (!fmpz_equal(a, b))
        {
            flint_printf("Fail: sparse dot sparse\n");
            abort();
        }

        fmpz_sparse_vec_dot_dense(a, u, x);

        if (!fmpz_equal(a, b))
        {
            flint_printf("Fail: sparse dot dense\n");
            abort();
        }
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_sparse_vec_clear(u);
        fmpz_sparse_vec_clear(v);
        _fmpz_vec_clear(w, len);
        _fmpz_vec_clear(x, len);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
