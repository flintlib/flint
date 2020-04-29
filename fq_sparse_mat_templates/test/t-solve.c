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
/* #include <sys/time.h> */

int
main(void)
{
    int iter, ret;
    slong rep, nreps = 100, r, c, i;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A, At;
    TEMPLATE(T, struct) *x, *x2, *b, *Atb, *Ax, *AtAx;
    slong niters[6] = {0, 0, 0, 0, 0, 0};
    slong psol[6] = {0, 0, 0, 0, 0, 0};
    slong nosol[6] = {0, 0, 0, 0, 0, 0};
    /* double elapsed[6] = {0, 0, 0, 0, 0}; */
    char *names[6] = {"rref", "lu", "Lanczos", "block Lanczos", "Wiedemann", "block Wiedemann"};
    /* struct timeval start, end; */
    FLINT_TEST_INIT(state);
    
    flint_printf("solving Ax = b....");
    fflush(stdout);
    
    for (rep = 0; rep < nreps; rep++)
    {
        if (rep % 5==0) {flint_printf("."); fflush(stdout);}
        TEMPLATE(T, ctx_randtest) (ctx, state);

        c = r = 50 + n_randint(state, 10);

        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        x = _TEMPLATE(T, vec_init) (c, ctx);
        x2 = _TEMPLATE(T, vec_init) (c, ctx);
        b = _TEMPLATE(T, vec_init) (r, ctx);
        Ax = _TEMPLATE(T, vec_init) (r, ctx);

        TEMPLATE(T, sparse_mat_randtest) (A, state, c/10, c/5, ctx);
        _TEMPLATE(T, vec_randtest) (x, state, c, ctx);
        TEMPLATE(T, sparse_mat_mul_vec) (b, A, x, ctx);

        for (i = 0; i < 6; ++i)
        {
            iter = 0;
            /* gettimeofday(&start, NULL); */
            switch (i) 
            {
            case 0: ret = TEMPLATE(T, sparse_mat_solve_rref) (x2, A, b, ctx); break;
            case 1: ret = TEMPLATE(T, sparse_mat_solve_lu) (x2, A, b, ctx); break;
            case 2: do ret = TEMPLATE(T, sparse_mat_solve_lanczos) (x2, A, b, state, ctx); while (ret == 0 && ++iter < 3); break;
            case 3: do ret = TEMPLATE(T, sparse_mat_solve_block_lanczos) (x2, A, b, 8, state, ctx); while (ret == 0 && ++iter < 3); break;
            case 4: ret = TEMPLATE(T, sparse_mat_solve_wiedemann) (x2, A, b, ctx); break;
            case 5: do ret = TEMPLATE(T, sparse_mat_solve_block_wiedemann) (x2, A, b, 8, state, ctx); while (ret == 0 && ++iter < 3); break;
            }
            /* gettimeofday(&end, NULL);
            elapsed[i] += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec); */
            if (ret == 0) nosol[i] += 1;
            else 
            {
                niters[i] += iter;
                TEMPLATE(T, sparse_mat_mul_vec) (Ax, A, x2, ctx);
                if (!_TEMPLATE(T, vec_equal) (b, Ax, A->r, ctx)) 
                {
                    if (i == 2 || i == 3)
                    {
                        TEMPLATE(T, sparse_mat_init) (At, c, r, ctx);
                        AtAx = _TEMPLATE(T, vec_init) (c, ctx);
                        Atb = _TEMPLATE(T, vec_init) (c, ctx);
                        TEMPLATE(T, sparse_mat_transpose) (At, A, ctx);
                        TEMPLATE(T, sparse_mat_mul_vec) (AtAx, At, Ax, ctx);
                        TEMPLATE(T, sparse_mat_mul_vec) (Atb, At, b, ctx);
                        if (!_TEMPLATE(T, vec_equal) (AtAx, Atb, A->c, ctx))
                        {
                            flint_printf("FAIL on %s: AtAx != Atb\n", names[i]);
                            abort();
                        } 
                        else psol[i] += 1;
                        _TEMPLATE(T, vec_clear) (AtAx, c, ctx);
                        _TEMPLATE(T, vec_clear) (Atb, c, ctx);
                        TEMPLATE(T, sparse_mat_clear) (At, ctx);
                    }
                    else
                    {
                        flint_printf("FAIL on %s: Ax != b\n", names[i]);
                        abort();
                    }
                }
            }
        }
        _TEMPLATE(T, vec_clear) (x, c, ctx);
        _TEMPLATE(T, vec_clear) (x2, c, ctx);
        _TEMPLATE(T, vec_clear) (b, r, ctx);
        _TEMPLATE(T, vec_clear) (Ax, r, ctx);
        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    for (i = 0; i < 6; ++i)
    {
        flint_printf("Solved with %s\n", names[i]);
        /* flint_printf("\tAverage time: %lf\n", elapsed[i]/nreps); */
        if (nosol[i])
            flint_printf("\tFound no solution for %wd/%wd examples\n", nosol[i], nreps);
        if (psol[i])    
            flint_printf("\tFound pseudo-solution for %wd/%wd examples\n", psol[i], nreps);
        if (niters[i])
            flint_printf("\tRequired %f extra iters per solution (on average).\n", (double)niters[i]/nreps);
    }
    return 0;
}

#endif
