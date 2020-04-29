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

slong TEMPLATE(T, sparse_mat_nullspace_rref) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, rk, numc, *Q;
    TEMPLATE(T, struct) *Xrow;
    TEMPLATE(T, sparse_mat_t) R;
    TEMPLATE(T, sparse_vec_struct) *Rrow;
    TEMPLATE(T, sparse_mat_init) (R, M->r, M->c, ctx);
    TEMPLATE(T, sparse_mat_set) (R, M, ctx);
    rk = TEMPLATE(T, sparse_mat_rref) (R, ctx);
    TEMPLATE(T, mat_init) (X, M->c, M->c-rk, ctx);
    if (rk != M->c) 
    {
        numc = 0;
        /* Mark which cols are pivots and enumerate the nonpivots */
        Q = flint_calloc(M->c, sizeof(*Q));
        for (i = 0; i < rk; ++i) 
            Q[R->rows[i].entries->ind] = -1;
        for (i = 0; i < M->c; ++i)
            if (Q[i]==UWORD(0)) Q[i] = numc++, TEMPLATE(T, one) (&X->rows[i][Q[i]], ctx);

        /* For each pivot col, set the corresponding row in X as */
        /* the negative of the associated row in R (reordered by Q) */
        for (i = 0; i < rk; ++i) 
        {
            Rrow = &R->rows[i];
            Xrow = X->rows[Rrow->entries[0].ind];
            for (j = 1; j < Rrow->nnz; ++j) 
                TEMPLATE(T, neg) (&Xrow[Q[Rrow->entries[j].ind]], Rrow->entries[j].val, ctx);
        }
        flint_free(Q);
    }
    TEMPLATE(T, sparse_mat_clear) (R, ctx);
    return X->c;
}

#endif
