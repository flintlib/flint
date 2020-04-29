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
    nmod_sparse_mat_t A, Ai;
    nmod_mat_t dA, dAiA;
    FLINT_TEST_INIT(state);
    
    flint_printf("inverting A....");
    fflush(stdout);
    
    for (rep = 0; rep < 100; rep++)
    {
        if (rep % 5==0) {flint_printf("."); fflush(stdout);}
        do r = n_randint(state, 200), c = n_randint(state, 200);
        while (r==UWORD(0) || c==UWORD(0));
        
        do n = n_randtest_not_zero(state);
        while (n <= 32 || !n_is_prime(n));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_sparse_mat_init(Ai, r, r, mod);
        nmod_sparse_mat_randtest(A, state, 0, c);
        nmod_mat_init(dA, r, c, n);
        nmod_mat_init(dAiA, r, c, n);
        nmod_sparse_mat_to_dense(dA, A);

        nmod_sparse_mat_inv(Ai, A);
        nmod_sparse_mat_mul_mat(dAiA, Ai, dA);
        nmod_mat_rref(dA);
        if (!nmod_mat_equal(dAiA, dA)) 
        {
            flint_printf("FAIL!\n");
            flint_printf("A^-1 x A = ");
            nmod_mat_print_pretty(dAiA);
            flint_printf("rref(A) = ");
            nmod_mat_print_pretty(dA);
            abort();
        }

        nmod_sparse_mat_clear(Ai);
        nmod_mat_clear(dA);
        nmod_mat_clear(dAiA);
        nmod_sparse_mat_clear(A);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
