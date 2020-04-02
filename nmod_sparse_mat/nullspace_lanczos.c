/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"
#include "nmod_sparse_mat.h"

slong nmod_sparse_mat_nullspace_lanczos(nmod_mat_t X, const nmod_sparse_mat_t A, flint_rand_t state, slong max_iters)
{
    /* Generate random solutions to a random system Ax=b and stop when nullspace filled */
    int ret;
    slong i, j, iter, nxs, *xps;
    mp_ptr x, x2, b, Ax, *xs;
    x = _nmod_vec_init(A->c);
    x2 = _nmod_vec_init(A->c);
    b = _nmod_vec_init(A->r);
    Ax = _nmod_vec_init(A->r);
    nxs = 0;
    xs = NULL;
    xps = NULL;
    for(iter = 0; iter < max_iters; )
    {
        _nmod_vec_randtest(x, state, A->c, A->mod);
        nmod_sparse_mat_mul_vec(b, A, x);
        if(nmod_sparse_mat_solve_lanczos(x2, A, b, state) == 0) {++iter; continue;} /* Lanczos failed */
        _nmod_vec_sub(x, x, x2, A->c, A->mod);
        if(_nmod_vec_is_zero(x, A->c)) {++iter; continue;}; /* Probably no nullspace */
        nmod_sparse_mat_mul_vec(Ax, A, x);
        if(!_nmod_vec_is_zero(Ax, A->c)) {++iter; continue;} /* In preimage of nullspace of A^t */ 
        
        /* Reduce by existing kernel vectors */
        for (j = nxs-1; j >= 0; --j) {
            _nmod_vec_scalar_addmul_nmod(x, xs[j], A->c, nmod_neg(x[xps[j]], A->mod), A->mod);
        }

        /* Normalize last nonzero entry to 1 */
        for (i = A->c-1; i >= 0 && x[i] == UWORD(0); --i);
        if (i == -1) {++iter; continue;} /* x in span of xs, nullspace probably complete */
        _nmod_vec_scalar_mul_nmod(x, x, A->c, nmod_inv(x[i], A->mod), A->mod);

        /* Reduce previous vectors by this one */
        for (j = 0; j < nxs; ++j) {
            _nmod_vec_scalar_addmul_nmod(xs[j], x, A->c, nmod_neg(xs[j][i], A->mod), A->mod);
        }

        /* Insert into list of vectors in nullspace */
        xs = realloc(xs, (nxs+1)*sizeof(*xs));
        xps = realloc(xps, (nxs+1)*sizeof(*xps));
        for (j = 0; j < nxs && i > xps[j]; ++j);
        memmove(xs + j + 1, xs + j, (nxs - j)*sizeof(*xs));
        memmove(xps + j + 1, xps + j, (nxs - j)*sizeof(*xps));
        xps[j] = i;
        xs[j] = x;
        nxs += 1;
        x = _nmod_vec_init(A->c); /* New vector for next iteration */
        iter = 0;
    }
    flint_free(xps);
    flint_free(x);
    flint_free(x2);
    flint_free(b);
    flint_free(Ax);
    nmod_mat_init(X, A->c, nxs, A->mod.n);
    for (i = 0; i < nxs; ++i)
    {
        for (j = 0; j < A->c; ++j) 
            X->rows[j][i] = xs[i][j];
        flint_free(xs[i]);
    }
    flint_free(xs);
    return X->c;
}
