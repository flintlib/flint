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
    slong rep, r, c;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_mat_t A, B;
    nmod_mat_t C, D;
    FLINT_TEST_INIT(state);
    

    flint_printf("conversion to/from dense matrix....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        r = n_randint(state, 10);
        c = n_randint(state, 10);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_sparse_mat_init(B, r, c, mod);
        nmod_mat_init(C, r, c, n);
        nmod_mat_init(D, r, c, n);
        
        nmod_sparse_mat_randtest(A, state, 0, c);
        nmod_sparse_mat_to_dense(C, A);
        nmod_sparse_mat_from_dense(B, C);
        
        if (!nmod_sparse_mat_equal(A, B))
        {
            flint_printf("FAIL: A != B\n");
            abort();
        }

        nmod_mat_randtest(C, state);
        nmod_sparse_mat_from_dense(A, C);
        nmod_sparse_mat_to_dense(D, A);
        
        if (!nmod_mat_equal(C, D))
        {
            flint_printf("FAIL: C != D\n");
            abort();
        }
        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
