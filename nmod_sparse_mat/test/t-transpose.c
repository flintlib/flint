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
    slong m, n, mod, mod2, rep;
    FLINT_TEST_INIT(state);
    

    flint_printf("transpose....");
    fflush(stdout);

    /* Rectangular transpose, same modulus */
    for (rep = 0; rep < 1000; rep++)
    {
        nmod_sparse_mat_t A, B, C;

        m = n_randint(state, 40);

        do
        mod = n_randtest_not_zero(state);
        while(mod <= 1);

        nmod_sparse_mat_init(A, m, mod);
        nmod_sparse_mat_randtest(A, state);
        
        nmod_sparse_mat_init(B, A->c, mod);
        nmod_sparse_mat_transpose(B, A);
        
        nmod_sparse_mat_init(C, m, mod);
        nmod_sparse_mat_transpose(C, B);
        
        if (!nmod_sparse_mat_equal(C, A))
        {
            flint_printf("FAIL: C != A\n");
            abort();
        }

        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(B);
        nmod_sparse_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
