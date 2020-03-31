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
    int niters = 0, nosol = 0, psolved = 0;
    slong rep, r, c, i;
    mp_limb_t n, a;
    nmod_t mod;
    nmod_sparse_mat_t A, At;
    mp_ptr x, x2, b, Atb, Ax, AtAx;
    double l_elapsed = 0, d_elapsed = 0;
    struct timeval start, end;
    FLINT_TEST_INIT(state);
    
    flint_printf("solving Ax=b....");
    fflush(stdout);
    
    for (rep = 0; rep < 100; rep++)
    {
        if(rep % 5==0) {flint_printf("."); fflush(stdout);}
        r = n_randint(state, 200);
        c = n_randint(state, 200);
        do n = n_randtest_not_zero(state);
        while (n <= 32 || !n_is_prime(n));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_sparse_mat_init(At, c, r, mod);

        nmod_sparse_mat_randtest(A, state, 0, c);
        nmod_sparse_mat_transpose(At, A);
        x = _nmod_vec_init(c);
        x2 = _nmod_vec_init(c);
        b = _nmod_vec_init(r);
        Ax = _nmod_vec_init(r);
        AtAx = _nmod_vec_init(c);
        Atb = _nmod_vec_init(c);

        _nmod_vec_randtest(x, state, c, mod);
        nmod_sparse_mat_mul_vec(b, A, x);

        /* Solve directly */
        gettimeofday(&start, NULL);
        ret = nmod_sparse_mat_solve_lu(x2, A, b);
        gettimeofday(&end, NULL);
        d_elapsed += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);
        nmod_sparse_mat_mul_vec(Ax, A, x2);
        if(!_nmod_vec_equal(b, Ax, A->r))
        {
            flint_printf("FAIL: Ax != b, got ret %d\n", ret);
            abort();
        } 

        /* Solve iteratively */
        gettimeofday(&start, NULL);
        for (iter=1; iter<=10; ++iter) /* TODO: set number of trials based on params */
            if (ret=nmod_sparse_mat_solve_lanczos(x2, A, b, state)) break;
        gettimeofday(&end, NULL);
        l_elapsed += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);

        if (iter==11)
        {
            nosol += 1;
            continue;
        }
        niters += iter;
        nmod_sparse_mat_mul_vec(Ax, A, x2);
        nmod_sparse_mat_mul_vec(AtAx, At, Ax);
        nmod_sparse_mat_mul_vec(Atb, At, b);
        if (!_nmod_vec_equal(AtAx, Atb, A->c))
        {
            flint_printf("FAIL: AtAx != Atb, got ret %d\n", ret);
            abort();
        } 
        else if (!_nmod_vec_equal(b, Ax, A->r))
        {
            psolved += 1;
        } 

        flint_free(x);
        flint_free(x2);
        flint_free(b);
        flint_free(Ax);
        flint_free(AtAx);
        flint_free(Atb);
        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(At);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    flint_printf("Average time for Lanzcos: %lf\n", l_elapsed/1000);
    flint_printf("Average time for direct: %lf\n", d_elapsed/1000);
    flint_printf("Lanczos found no solution for %wd/%wd examples, pseudo-solution for %wd/%wd examples, and required %f iters per solution.\n", nosol, 1000, psolved, 1000, niters/1000.);
    return 0;
}
