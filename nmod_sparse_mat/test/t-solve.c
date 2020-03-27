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
    slong m, n, mod, rep;
    FLINT_TEST_INIT(state);
    

    flint_printf("solving Ax=b....");
    fflush(stdout);
    int niters = 0, nosol = 0, psolved = 0, nusolved = 0;

    for (rep = 0; rep < 1000; rep++)
    {
        nmod_sparse_mat_t A, At;
        mp_ptr x, x2, b, Atb, Ax, AtAx;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        do
            mod = n_randtest_not_zero(state);
        while(mod <= 32 || !n_is_prime(mod));

        nmod_sparse_mat_init(A, m, mod);
        nmod_sparse_mat_randtest(A, state);
        nmod_sparse_mat_init(At, A->c, mod);
        nmod_sparse_mat_transpose(At, A);
        x = _nmod_vec_init(A->c);
        x2 = _nmod_vec_init(A->c);
        b = _nmod_vec_init(A->r);
        Ax = _nmod_vec_init(A->r);
        AtAx = _nmod_vec_init(A->c);
        Atb = _nmod_vec_init(A->c);

        _nmod_vec_randtest(x, state, A->c, A->mod);
        nmod_sparse_mat_mul_vec(b, A, x);
        nmod_sparse_mat_mul_vec(Atb, At, b);
        int iter, ret;
        for(iter=1; iter<=10; ++iter) 
            if(ret=nmod_sparse_mat_solve_lanczos(x2, A, b, state)) break;
        if(iter==11) {
            nosol += 1;
            continue;
        }
        niters += iter;
        nmod_sparse_mat_mul_vec(Ax, A, x2);
        nmod_sparse_mat_mul_vec(AtAx, At, Ax);
        if (!_nmod_vec_equal(AtAx, Atb, A->c))
        {
            flint_printf("FAIL: AtAx != Atb for mod=%wd, got ret %d\n", mod, ret);
            abort();
        } else if(!_nmod_vec_equal(b, Ax, A->r)) {
            psolved += 1;
        } else if(!_nmod_vec_equal(x, x2, A->c)) {
            nusolved += 1;
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
    flint_printf("No solution found for %wd/%wd examples\n", nosol, 1000);
    flint_printf("Average number of iters to find solution: %f\n", niters/1000.);
    flint_printf("Pseudo-solution found for %wd/%wd examples\n", psolved, 1000);
    flint_printf("Alternate solution found for %wd/%wd examples\n", nusolved, 1000);
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
