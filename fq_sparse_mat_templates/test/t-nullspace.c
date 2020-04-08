/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "ulong_extras.h"

int
main(void)
{
    int iter, ret;
    slong rep, nreps = 100, r, c, i;
    TEMPLATE(T, t) a;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A;
    TEMPLATE(T, mat_t) X, AX;
    slong rk[4];
    slong discrep[4] = {0, 0, 0, 0};
    double elapsed[4] = {0, 0, 0, 0};
    struct timeval start, end;
    FLINT_TEST_INIT(state);
    
    flint_printf("finding nullspace of A....");
    fflush(stdout);
    
    for (rep = 0; rep < nreps; )
    {
        if (rep % 5==0) {flint_printf("."); fflush(stdout);}
        TEMPLATE(T, ctx_randtest) (ctx, state);
        r = c = n_randint(state, 200);
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, sparse_mat_randtest) (A, state, c/5, c/2, ctx);
        for (i = 0; i < 4; ++i) 
        {
            gettimeofday(&start, NULL);
            switch (i) 
            {
            case 0: rk[0] = TEMPLATE(T, sparse_mat_nullspace_rref) (X, A, ctx); break;
            case 1: rk[1] = TEMPLATE(T, sparse_mat_nullspace_lu) (X, A, ctx); break;
            case 2: rk[2] = TEMPLATE(T, sparse_mat_nullspace_lanczos) (X, A, state, 5, ctx); break;
            case 3: rk[3] = TEMPLATE(T, sparse_mat_nullspace_wiedemann) (X, A, state, 5, ctx); break;
            }
            if (i == 0 && rk[0] == UWORD(0) ) { TEMPLATE(T, mat_clear) (X, ctx); break;}
            else ++rep;
            gettimeofday(&end, NULL);
            elapsed[i] += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);
            if (rk[i]!=0) 
            {
                TEMPLATE(T, mat_init) (AX, A->r, X->c, ctx);
                TEMPLATE(T, sparse_mat_mul_mat) (AX, A, X, ctx); 
                if (!TEMPLATE(T, mat_is_zero) (AX, ctx)) 
                {
                    flint_printf("FAIL: %d!\n", i);
                    if (i != 0)
                        flint_printf("Nullity should be %wd\n", rk[0]);
                    TEMPLATE(T, sparse_mat_print_pretty) (A, ctx);
                    TEMPLATE(T, mat_print_pretty) (X, ctx);
                    abort();
                }
                TEMPLATE(T, mat_clear) (AX, ctx);
            }
            if (rk[i] != rk[0]) 
            {
                discrep[i] += 1;
            }
            TEMPLATE(T, mat_clear) (X, ctx);
        }
        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    flint_printf("Average time for rref: %lf\n", elapsed[0]/nreps);
    flint_printf("Average time for LU: %lf\n", elapsed[1]/nreps);
    flint_printf("Average time for Lanzcos: %lf\n", elapsed[2]/nreps);
    flint_printf("Average time for Wiedemann: %lf\n", elapsed[3]/nreps);
    flint_printf("LU did not find full nullspace for %wd/%wd examples\n", discrep[1], nreps);
    flint_printf("Lanczos did not find full nullspace for %wd/%wd examples\n", discrep[2], nreps);
    flint_printf("Wiedemann did not find full nullspace for %wd/%wd examples\n", discrep[3], nreps);
    return 0;
}

#endif
