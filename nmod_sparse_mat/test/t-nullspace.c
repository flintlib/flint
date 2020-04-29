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
/* #include <sys/time.h> */
#include "flint.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong rep, nreps = 100, r, c, i;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_mat_t A;
    nmod_mat_t X, AX;
    slong rk[6];
/*     double elapsed[6] = {0, 0, 0, 0, 0, 0}; */
    slong discrep[6] = {0, 0, 0, 0, 0, 0};
    char *names[6] = {"rref", "lu", "Lanczos", "block Lanczos", "Wiedemann", "block Wiedemann"};
/*     struct timeval start, end; */
    FLINT_TEST_INIT(state);
    
    flint_printf("finding nullspace of A....");
    fflush(stdout);
    
    for (rep = 0; rep < nreps; ++rep)
    {
        if (rep % 5==0) {flint_printf("."); fflush(stdout);}
        c = r = 200 + n_randint(state, 100);
        do n = n_randtest_not_zero(state);
        while (n <= 32 || !n_is_prime(n));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_sparse_mat_randtest(A, state, c/20, c/10);
        for (i = 0; i < 6; ++i) 
        {
/*             gettimeofday(&start, NULL); */
            switch (i) 
            {
            case 0: rk[0] = nmod_sparse_mat_nullspace_rref(X, A); break;
            case 1: rk[1] = nmod_sparse_mat_nullspace_lu(X, A); break;
            case 2: rk[2] = nmod_sparse_mat_nullspace_lanczos(X, A, state, 2); break;
            case 3: rk[3] = nmod_sparse_mat_nullspace_block_lanczos(X, A, 8, state, 2); break;
            case 4: rk[4] = nmod_sparse_mat_nullspace_wiedemann(X, A, state, 2); break;
            case 5: rk[5] = nmod_sparse_mat_nullspace_block_wiedemann(X, A, 8, state, 2); break;
            }
            /* if (i == 0 && rk[0] == UWORD(0) ) { nmod_mat_clear(X); break;} */
/*             gettimeofday(&end, NULL);
            elapsed[i] += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec); */
            if (X->c==0) continue;
            nmod_mat_init(AX, A->r, X->c, n);
            nmod_sparse_mat_mul_mat(AX, A, X); 
            /*nmod_mat_print_pretty(AX); */
            if (!nmod_mat_is_zero(AX)) 
            {
                flint_printf("FAIL: %s!\n", names[i]);
                if (i != 0)
                    flint_printf("Nullity should be %wd\n", rk[i-1]);
                nmod_sparse_mat_print_pretty(A);
                nmod_mat_print_pretty(X);
                abort();
            }
            if (rk[i] != rk[0]) discrep[i] += 1;
            
            nmod_mat_clear(X);
            nmod_mat_clear(AX);
        }
        nmod_sparse_mat_clear(A);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    for (i = 0; i < 6; ++i)
    {
        flint_printf("Found nullspace with %s\n", names[i]);
/*         flint_printf("\tAverage time %lf:\n", elapsed[i]/nreps); */
        if (discrep[i] > 0)
            flint_printf("\tFailed to find full nullspace in %wd/%wd trials\n", 
                         discrep[i], nreps);        
    }
    return 0;
}
