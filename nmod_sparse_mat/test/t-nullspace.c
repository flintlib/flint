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
    int iter, ret;
    slong rep, nreps = 100, r, c, i, lanczos_discrep = 0, wied_discrep = 0;
    mp_limb_t n, a;
    nmod_t mod;
    nmod_sparse_mat_t A;
    nmod_mat_t X, AX;
    slong rk[3];
    double elapsed[4] = {0, 0, 0, 0};
    struct timeval start, end;
    FLINT_TEST_INIT(state);
    
    flint_printf("finding nullspace of A....");
    fflush(stdout);
    
    for (rep = 0; rep < nreps; rep++)
    {
        if(rep % 5==0) {flint_printf("."); fflush(stdout);}
        r = n_randint(state, 500);
        c = r + n_randint(state, 10);
        do n = n_randtest_not_zero(state);
        while (n <= 32 || !n_is_prime(n));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_sparse_mat_randtest(A, state, c/20, c/10);
        for (i = 0; i < 4; ++i) 
        {
            gettimeofday(&start, NULL);
            switch (i) 
            {
            case 0: rk[0] = nmod_sparse_mat_nullspace_rref(X, A); break;
            case 1: rk[1] = nmod_sparse_mat_nullspace_lu(X, A); break;
            case 2: rk[2] = nmod_sparse_mat_nullspace_lanczos(X, A, state, 2); break;
            case 3: rk[3] = nmod_sparse_mat_nullspace_wiedemann(X, A, state, 2); break;
            }
            gettimeofday(&end, NULL);
            elapsed[i] += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);
            if(X->c==0) continue;
            nmod_mat_init(AX, A->r, X->c, n);
            nmod_sparse_mat_mul_mat(AX, A, X); 
            /*nmod_mat_print_pretty(AX); */
            if (!nmod_mat_is_zero(AX)) 
            {
                flint_printf("FAIL: %d!\n", i);
                if (i != 0)
                    flint_printf("Nullity should be %wd\n", rk[i-1]);
                nmod_sparse_mat_print_pretty(A);
                nmod_mat_print_pretty(X);
                abort();
            }
            if (i == 2 && rk[i] != rk[0])
            {
                lanczos_discrep += 1;
            }
            if (i == 3 && rk[i] != rk[0])
            {
                wied_discrep += 1;
            }
            
            nmod_mat_clear(X);
            nmod_mat_clear(AX);
        }
        nmod_sparse_mat_clear(A);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    flint_printf("Average time for rref: %lf\n", elapsed[0]/nreps);
    flint_printf("Average time for LU: %lf\n", elapsed[1]/nreps);
    flint_printf("Average time for Lanzcos: %lf\n", elapsed[2]/nreps);
    flint_printf("Average time for Wiedemann: %lf\n", elapsed[3]/nreps);
    flint_printf("Lanczos did not find full nullspace for %wd/%wd examples\n", lanczos_discrep, nreps);
    flint_printf("Wiedemann did not find full nullspace for %wd/%wd examples\n", lanczos_discrep, nreps);
    return 0;
}
