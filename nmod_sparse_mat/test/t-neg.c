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
    

    flint_printf("neg....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        nmod_sparse_mat_t A, B, C, D;

        m = n_randint(state, 200);
        do
            mod = n_randlimb(state);
        while(mod <= UWORD(1));

        nmod_sparse_mat_init(A, m, mod);
        nmod_sparse_mat_init(B, m, mod);
        nmod_sparse_mat_init(C, m, mod);
        nmod_sparse_mat_init(D, m, mod);

        nmod_sparse_mat_randtest(A, state);
        nmod_sparse_mat_randtest(B, state);

        nmod_sparse_mat_sub(C, A, B);
        nmod_sparse_mat_neg(B, B);
        nmod_sparse_mat_add(D, A, B);

        if (!nmod_sparse_mat_equal(C, D))
        {
            flint_printf("FAIL\n");
            abort();
        }

        nmod_sparse_mat_neg(C, B);
        nmod_sparse_mat_neg(B, B);

        if (!nmod_sparse_mat_equal(C, B))
        {
            flint_printf("FAIL\n");
            abort();
        }
        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(B);
        nmod_sparse_mat_clear(C);
        nmod_sparse_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
