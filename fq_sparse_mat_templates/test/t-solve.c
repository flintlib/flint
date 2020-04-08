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
#include "ulong_extras.h"
#include <sys/time.h>

int
main(void)
{
    int iter, ret;
    int niters = 0, wied_nosol = 0, nosol = 0, psolved = 0;
    slong rep, r, c, i, nrep = 100;
    TEMPLATE(T, t) a;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A, At;
    TEMPLATE(T, struct) *x, *x2, *b, *Atb, *Ax, *AtAx;
    double rref_elapsed = 0, lu_elapsed = 0, lanczos_elapsed = 0, wiedemann_elapsed;
    struct timeval start, end;
    FLINT_TEST_INIT(state);
    
    flint_printf("solving Ax=b....");
    fflush(stdout);
    
    for (rep = 0; rep < nrep; rep++)
    {
        if(rep % 5==0) {flint_printf("."); fflush(stdout);}
        TEMPLATE(T, ctx_randtest) (ctx, state);

        do c = r = n_randint(state, 200);
        while(c == 0 || r == 0);

        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (At, c, r, ctx);
        x = _TEMPLATE(T, vec_init) (c, ctx);
        x2 = _TEMPLATE(T, vec_init) (c, ctx);
        b = _TEMPLATE(T, vec_init) (r, ctx);
        Ax = _TEMPLATE(T, vec_init) (r, ctx);
        AtAx = _TEMPLATE(T, vec_init) (c, ctx);
        Atb = _TEMPLATE(T, vec_init) (c, ctx);

        TEMPLATE(T, sparse_mat_randtest) (A, state, c/20, c/10, ctx);
        TEMPLATE(T, sparse_mat_transpose) (At, A, ctx);

        _TEMPLATE(T, vec_randtest) (x, state, c, ctx);
        TEMPLATE(T, sparse_mat_mul_vec) (b, A, x, ctx);
        /* Solve via reduced row echelon form */
         gettimeofday(&start, NULL);
        ret = TEMPLATE(T, sparse_mat_solve_rref) (x2, A, b, ctx);
        gettimeofday(&end, NULL);
        rref_elapsed += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);
        TEMPLATE(T, sparse_mat_mul_vec) (Ax, A, x2, ctx);
        if(!_TEMPLATE(T, vec_equal) (b, Ax, A->r, ctx))
        {
            flint_printf("FAIL: Ax != b, got ret %d\n", ret);
            abort();
        } 

        /* Solve via lu decomposition */
        gettimeofday(&start, NULL);
        ret = TEMPLATE(T, sparse_mat_solve_lu) (x2, A, b, ctx);
        gettimeofday(&end, NULL);
        lu_elapsed += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);
        TEMPLATE(T, sparse_mat_mul_vec) (Ax, A, x2, ctx);
        if(!_TEMPLATE(T, vec_equal) (b, Ax, A->r, ctx))
        {
            flint_printf("FAIL: Ax != b, got ret %d\n", ret);
            abort();
        } 

        /* Solve iteratively */
/**/         gettimeofday(&start, NULL);
        ret=TEMPLATE(T, sparse_mat_solve_wiedemann) (x2, A, b, ctx);
        gettimeofday(&end, NULL);
        wiedemann_elapsed += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);
        if (ret == 0)
        {
            wied_nosol += 1;
        }
        else
        {
            TEMPLATE(T, sparse_mat_mul_vec) (Ax, A, x2, ctx);
            if (!_TEMPLATE(T, vec_equal) (b, Ax, A->r, ctx))
            {
                flint_printf("FAIL: Ax != b\n");
                abort();
            }
        }

        gettimeofday(&start, NULL);
        iter = 0;
        do ret=TEMPLATE(T, sparse_mat_solve_lanczos) (x2, A, b, state, ctx);
        while(ret==0 && ++iter < 30);
        gettimeofday(&end, NULL);
        lanczos_elapsed += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec);
        if (ret==0)
        {
            nosol += 1;
            continue;
        }
        else
        {
            niters += iter;
            TEMPLATE(T, sparse_mat_mul_vec) (Ax, A, x2, ctx);
            TEMPLATE(T, sparse_mat_mul_vec) (AtAx, At, Ax, ctx);
            TEMPLATE(T, sparse_mat_mul_vec) (Atb, At, b, ctx);
            if (!_TEMPLATE(T, vec_equal) (AtAx, Atb, A->c, ctx))
            {
                flint_printf("FAIL: AtAx != Atb, got ret %d\n", ret);
                abort();
            } 
            else if (!_TEMPLATE(T, vec_equal) (b, Ax, A->r, ctx))
            {
                psolved += 1;
            }
        }
 
        _TEMPLATE(T, vec_clear) (x, c, ctx);
        _TEMPLATE(T, vec_clear) (x2, c, ctx);
        _TEMPLATE(T, vec_clear) (b, r, ctx);
        _TEMPLATE(T, vec_clear) (Ax, r, ctx);
        _TEMPLATE(T, vec_clear) (AtAx, c, ctx);
        _TEMPLATE(T, vec_clear) (Atb, c, ctx);
        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, sparse_mat_clear) (At, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    flint_printf("Average time for Wiedemann: %lf\n", wiedemann_elapsed/nrep);
    flint_printf("Average time for Lanzcos: %lf\n", lanczos_elapsed/nrep);
    flint_printf("Average time for LU: %lf\n", lu_elapsed/nrep);
    flint_printf("Average time for rref: %lf\n", rref_elapsed/nrep);
    flint_printf("Wiedemann found no solution for %wd/%wd examples.\n", wied_nosol, nrep);
    flint_printf("Lanczos found no solution for %wd/%wd examples, pseudo-solution for %wd/%wd examples, and required %f extra iters per solution (on average).\n", nosol, nrep, psolved, nrep, (double)niters/nrep);
    return 0;
}

#endif
