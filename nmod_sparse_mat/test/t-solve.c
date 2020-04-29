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
/*#include <sys/time.h>*/
#include "flint.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int iter, ret;
    slong rep, nreps = 100, r, c, i;
    mp_limb_t n;
    nmod_t mod;
    nmod_mat_t dA;
    nmod_sparse_mat_t A, At;
    mp_ptr x, x2, b, Atb, Ax, AtAx;
    slong niters[6] = {0, 0, 0, 0, 0, 0};
    slong psol[6] = {0, 0, 0, 0, 0, 0};
    slong nosol[6] = {0, 0, 0, 0, 0, 0};
    /*double elapsed[6] = {0, 0, 0, 0, 0};*/
    char *names[6] = {"rref", "lu", "Lanczos", "block Lanczos", "Wiedemann", "block Wiedemann"};
    /*struct timeval start, end;*/
    FLINT_TEST_INIT(state);
    
    flint_printf("solving Ax = b....");
    fflush(stdout);
    
    for (rep = 0; rep < nreps; rep++)
    {
        if (rep % 5==0) {flint_printf("."); fflush(stdout);}
        
        c = r = 100 + n_randint(state, 50);
        do n = n_randtest_not_zero(state);
        while (n <= 32 || !n_is_prime(n));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_mat_init(dA, r, c, n);
        x = _nmod_vec_init(c);
        x2 = _nmod_vec_init(c);
        b = _nmod_vec_init(r);
        Ax = _nmod_vec_init(r);

        nmod_sparse_mat_randtest(A, state, c/20, c/10);
        nmod_sparse_mat_to_dense(dA, A);
        _nmod_vec_randtest(x, state, c, mod);
        nmod_sparse_mat_mul_vec(b, A, x);

        for (i = 0; i < 6; ++i)
        {
            iter = 0;
            /*gettimeofday(&start, NULL);*/
            switch (i) 
            {
            case 0: ret = nmod_sparse_mat_solve_rref(x2, A, b); break;
            case 1: ret = nmod_sparse_mat_solve_lu(x2, A, b); break;
            case 2: do ret = nmod_sparse_mat_solve_lanczos(x2, A, b, state); while (ret == 0 && ++iter < 30); break;
            case 3: do ret = nmod_sparse_mat_solve_block_lanczos(x2, A, b, 8, state); while (ret == 0 && ++iter < 30); break;
            case 4: ret = nmod_sparse_mat_solve_wiedemann(x2, A, b); break;
            case 5: do ret = nmod_sparse_mat_solve_block_wiedemann(x2, A, b, 4, state); while (ret == 0 && ++iter < 3); break;
            }
            /*gettimeofday(&end, NULL);
            elapsed[i] += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);*/
            if (ret == 0) nosol[i] += 1;
            else 
            {
                niters[i] += iter;
                nmod_sparse_mat_mul_vec(Ax, A, x2);
                if (!_nmod_vec_equal(b, Ax, A->r)) 
                {
                    if (i == 2 || i == 3)
                    {
                        nmod_sparse_mat_init(At, c, r, mod);
                        AtAx = _nmod_vec_init(c);
                        Atb = _nmod_vec_init(c);
                        nmod_sparse_mat_transpose(At, A);
                        nmod_sparse_mat_mul_vec(AtAx, At, Ax);
                        nmod_sparse_mat_mul_vec(Atb, At, b);
                        if (!_nmod_vec_equal(AtAx, Atb, A->c))
                        {
                            flint_printf("FAIL on %s: AtAx != Atb\n", names[i]);
                            abort();
                        } 
                        else psol[i] += 1;
                        flint_free(AtAx);
                        flint_free(Atb);
                        nmod_sparse_mat_clear(At);
                    }
                    else
                    {
                        flint_printf("FAIL on %s: Ax != b\n", names[i]);
                        flint_printf("A = ");
                        nmod_sparse_mat_print_pretty(A);
                        flint_printf("x = ");
                        _nmod_vec_print_pretty(x, c, mod);
                        flint_printf("x2 = ");
                        _nmod_vec_print_pretty(x2, c, mod);
                        flint_printf("Ax2 = ");
                        _nmod_vec_print_pretty(Ax, r, mod);
                        flint_printf("b = ");
                        _nmod_vec_print_pretty(b, r, mod);
                        abort();
                    }
                }
            }
        }
        flint_free(x);
        flint_free(x2);
        flint_free(b);
        flint_free(Ax);
        nmod_sparse_mat_clear(A);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    for (i = 0; i < 6; ++i)
    {
        flint_printf("Solved with %s\n", names[i]);
        /*flint_printf("\tAverage time: %lf\n", elapsed[i]/nreps);*/
        if (nosol[i])
            flint_printf("\tFound no solution for %wd/%wd examples\n", nosol[i], nreps);
        if (psol[i])    
            flint_printf("\tFound pseudo-solution for %wd/%wd examples\n", psol[i], nreps);
        if (niters[i])
            flint_printf("\tRequired %f extra iters per solution (on average).\n", (double)niters[i]/nreps);
    }
    return 0;
}
