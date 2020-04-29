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

slong nmod_sparse_mat_nullspace_block_wiedemann(nmod_mat_t X, const nmod_sparse_mat_t M, slong block_size, flint_rand_t state, slong max_iters)
{
    /* Generate random solutions to a random system Mx = b and stop when nullspace filled */
    slong i, j, iter, nxs, *xps;
    mp_ptr x, *xs;
    x = _nmod_vec_init(M->c);
    nxs = 0;
    xs = NULL;
    xps = NULL;
    for (iter = 0; iter < max_iters; )
    {
        if (nmod_sparse_mat_nullvector_block_wiedemann(x, M, block_size, state) == 0) {++iter; continue;}
        
        /* Reduce by existing kernel vectors */
        for (j = nxs-1; j >= 0; --j) 
        {
            _nmod_vec_scalar_addmul_nmod(x, xs[j], M->c, nmod_neg(x[xps[j]], M->mod), M->mod);
        }

        /* Normalize last nonzero entry to 1 */
        for (i = M->c-1; i >= 0 && x[i] == UWORD(0); --i);
        if (i == -1) {++iter; continue;} /* x in span of xs, nullspace probably complete */
        _nmod_vec_scalar_mul_nmod(x, x, M->c, nmod_inv(x[i], M->mod), M->mod);

        /* Reduce previous vectors by this one */
        for (j = 0; j < nxs; ++j) 
        {
            _nmod_vec_scalar_addmul_nmod(xs[j], x, M->c, nmod_neg(xs[j][i], M->mod), M->mod);
        }

        /* Insert into list of vectors in nullspace (ordered by pivot) */
        xs = realloc(xs, (nxs+1)*sizeof(*xs));
        xps = realloc(xps, (nxs+1)*sizeof(*xps));
        for (j = 0; j < nxs && i > xps[j]; ++j);
        memmove(xs + j + 1, xs + j, (nxs - j)*sizeof(*xs));
        memmove(xps + j + 1, xps + j, (nxs - j)*sizeof(*xps));
        xps[j] = i;
        xs[j] = x;
        nxs += 1;
        x = _nmod_vec_init(M->c); /* New vector for next iteration */
        iter = 0;
    }
    flint_free(xps);
    flint_free(x);
    nmod_mat_init(X, M->c, nxs, M->mod.n);
    for (i = 0; i < nxs; ++i)
    {
        for (j = 0; j < M->c; ++j) 
            X->rows[j][i] = xs[i][j];
        flint_free(xs[i]);
    }
    flint_free(xs);
    return X->c;
}
