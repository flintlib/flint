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
#include <sys/time.h>
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
    nmod_mat_t dA;
    FLINT_TEST_INIT(state);
    
    flint_printf("converting A to reduced row echelon form....");
    fflush(stdout);
    
    for (rep = 0; rep < 1000; rep++)
    {
        r = n_randint(state, 100);
        c = n_randint(state, 100);
        do n = n_randtest_not_zero(state);
        while (n <= 32 || !n_is_prime(n));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_sparse_mat_init(B, r, c, mod);
        nmod_mat_init(dA, r, c, n);

        nmod_sparse_mat_randtest(A, state, 0, n);
        nmod_sparse_mat_to_dense(dA, A);
        nmod_sparse_mat_rref(A);
        nmod_mat_rref(dA);
        nmod_sparse_mat_from_dense(B, dA);
        if (!nmod_sparse_mat_equal(A, B)) 
        {
            flint_printf("FAIL!\n");
            abort();
        }

        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(B);
        nmod_mat_clear(dA);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
