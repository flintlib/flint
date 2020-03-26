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
    slong m, mod, rep;
    FLINT_TEST_INIT(state);
    

    flint_printf("conversion to/from dense matrix....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        nmod_sparse_mat_t A, B;
        nmod_mat_t C, D;

        m = n_randint(state, 20);

        do
            mod = n_randtest_not_zero(state);
        while(mod <= 1);

        nmod_sparse_mat_init(A, m, mod);
        nmod_sparse_mat_randtest(A, state);

        
        nmod_mat_init(C, m, A->c, mod);
        nmod_sparse_mat_to_dense(C, A);

        nmod_sparse_mat_init(B, m, mod);
        nmod_sparse_mat_randtest(B, state);
        nmod_sparse_mat_from_dense(B, C);
        
        if (!nmod_sparse_mat_equal(A, B))
        {
            flint_printf("FAIL: A != B\n");
            abort();
        }

        nmod_mat_randtest(C, state);
        nmod_sparse_mat_from_dense(A, C);

        nmod_mat_init(D, m, C->c, mod);
        nmod_sparse_mat_to_dense(D, A);
        
        if(!nmod_mat_equal(C, D)) {
            flint_printf("FAIL: C != D\n");
            abort();
        }
        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(B);
        nmod_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
