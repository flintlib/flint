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
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, mod, rep;
    FLINT_TEST_INIT(state);
    

    flint_printf("multiplication by vec/mat....");
    fflush(stdout);

    for (rep = 0; rep < 10; rep++)
    {
        nmod_sparse_mat_t A;
        mp_ptr x, y, y2;
        nmod_mat_t B, X, Y, Y2;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        do
            mod = n_randtest_not_zero(state);
        while(mod <= 1);

        nmod_sparse_mat_init(A, m, mod);
        nmod_sparse_mat_randtest(A, state);
        nmod_mat_init(B, m, A->c, mod);
        nmod_sparse_mat_to_dense(B, A);

        x = _nmod_vec_init(A->c);
        y = _nmod_vec_init(A->r);
        y2 = _nmod_vec_init(A->r);
        _nmod_vec_randtest(x, state, A->c, A->mod);
        nmod_sparse_mat_mul_vec(y, A, x);
        nmod_mat_mul_vec(y2, B, x);

        if (!_nmod_vec_equal(y, y2, A->r))
        {
            flint_printf("FAIL: y != y2\n");
            abort();
        }

        nmod_mat_init(X, A->c, n, mod);
        nmod_mat_init(Y, A->r, n, mod);
        nmod_mat_init(Y2, A->r, n, mod);
        nmod_mat_randtest(X, state);
        nmod_sparse_mat_mul_mat(Y, A, X);
        nmod_mat_mul(Y2, B, X);

        if (!nmod_mat_equal(Y, Y2)) {
            flint_printf("Fail: Y != Y2\n");
            abort();
        }

        nmod_sparse_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(X);
        nmod_mat_clear(Y);
        nmod_mat_clear(Y2);
        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
        _nmod_vec_clear(y2);
     }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
