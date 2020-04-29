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
    slong rep, bits, r, c, k;
    fmpz_sparse_mat_t A;
    fmpz *x, *y, *y2;
    fmpz_mat_t B, X, Y, Y2;
    FLINT_TEST_INIT(state);
    

    flint_printf("multiplication by vec/mat....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 70);
        while (bits < UWORD(2));
        r = 1 + n_randint(state, 20);
        c = 1 + n_randint(state, 20);
        k = 1 + n_randint(state, 20);
        fmpz_sparse_mat_init(A, r, c);
        fmpz_mat_init(B, r, c);
        fmpz_mat_init(X, c, k);
        fmpz_mat_init(Y, r, k);
        fmpz_mat_init(Y2, r, k);
        x = _fmpz_vec_init(c);
        y = _fmpz_vec_init(r);
        y2 = _fmpz_vec_init(r);

        fmpz_sparse_mat_randtest(A, state, 0, c, bits);
        fmpz_sparse_mat_to_dense(B, A);
        _fmpz_vec_randtest(x, state, c, bits);
        fmpz_sparse_mat_mul_vec(y, A, x);
        fmpz_mat_mul_vec(y2, B, x);

        if (!_fmpz_vec_equal(y, y2, A->r))
        {
            flint_printf("FAIL: y != y2\n");
            abort();
        }

        fmpz_mat_randtest(X, state, bits);
        fmpz_sparse_mat_mul_mat(Y, A, X);
        fmpz_mat_mul(Y2, B, X);

        if (!fmpz_mat_equal(Y, Y2))
        {
            flint_printf("Fail: Y != Y2\n");
            abort();
        }

        fmpz_sparse_mat_clear(A);
        fmpz_mat_clear(B);
        _fmpz_vec_clear(x, c);
        _fmpz_vec_clear(y, r);
        _fmpz_vec_clear(y2, r);
        fmpz_mat_clear(X);
        fmpz_mat_clear(Y);
        fmpz_mat_clear(Y2);
     }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
