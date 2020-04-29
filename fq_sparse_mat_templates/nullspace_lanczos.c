/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include <string.h>
#include "templates.h"

slong TEMPLATE(T, sparse_mat_nullspace_lanczos) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, slong max_iters, const TEMPLATE(T, ctx_t) ctx)
{
    /* Generate random solutions to a random system Mx = b and stop when nullspace filled */
    slong i, j, iter, nxs, *xps;
    TEMPLATE(T, t) cc;
    TEMPLATE(T, struct) *x, **xs;

    TEMPLATE(T, init) (cc, ctx);
    x = _TEMPLATE(T, vec_init) (M->c, ctx);
    nxs = 0;
    xs = NULL;
    xps = NULL;
    for (iter = 0; iter < max_iters; )
    {
        if (TEMPLATE(T, sparse_mat_nullvector_lanczos) (x, M, state, ctx) == 0) {++iter; continue;}

        /* Reduce by existing kernel vectors */
        for (j = nxs-1; j >= 0; --j) 
        {
            TEMPLATE(T, neg) (cc, &x[xps[j]], ctx);
            _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (x, xs[j], M->c, cc, ctx);
        }

        /* Normalize last nonzero entry to 1 */
        for (i = M->c-1; i >= 0 && TEMPLATE(T, is_zero) (&x[i], ctx); --i);
        if (i == -1) {++iter; continue;} /* x in span of xs, nullspace probably complete */
        TEMPLATE(T, inv) (cc, &x[i], ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_mul, T)) (x, x, M->c, cc, ctx);

        /* Reduce previous vectors by this one */
        for (j = 0; j < nxs; ++j) 
        {
            TEMPLATE(T, neg) (cc, &xs[j][i], ctx);
            _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (xs[j], x, M->c, cc, ctx);
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
        x = _TEMPLATE(T, vec_init) (M->c, ctx); /* New vector for next iteration */
        iter = 0;
    }
    flint_free(xps);
    TEMPLATE(T, clear) (cc, ctx);
    _TEMPLATE(T, vec_clear) (x, M->c, ctx);
    TEMPLATE(T, mat_init) (X, M->c, nxs, ctx);
    for (i = 0; i < nxs; ++i)
    {
        for (j = 0; j < M->c; ++j) 
            TEMPLATE(T, set) (&X->rows[j][i], &xs[i][j], ctx);
        _TEMPLATE(T, vec_clear) (xs[i], M->c, ctx);
    }
    flint_free(xs);
    return X->c;
}

#endif
