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

slong TEMPLATE(T, sparse_mat_nullspace_lu) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx)
{
    slong rk, *P, *Q, *Qi, i, j;
    TEMPLATE(T, t) cc;
    TEMPLATE(T, sparse_mat_t) L, U;
    TEMPLATE(T, sparse_entry_struct) *e;
    TEMPLATE(T, sparse_vec_struct) *Urow;
    TEMPLATE(T, struct) *Xrow;

    P = flint_malloc(M->r * sizeof(*P));
    Q = flint_malloc(M->c * sizeof(*Q));
    TEMPLATE(T, init) (cc, ctx);
    TEMPLATE(T, sparse_mat_init) (L, M->r, M->c, ctx);
    TEMPLATE(T, sparse_mat_init) (U, M->r, M->c, ctx);
    rk = TEMPLATE(T, sparse_mat_lu) (P, Q, L, U, M, ctx);
    flint_free(P);
    TEMPLATE(T, sparse_mat_clear) (L, ctx);
    for (i = 0; i < rk; ++i)
    {
        TEMPLATE(T, inv) (cc, U->rows[i].entries[0].val, ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (&U->rows[i], &U->rows[i], cc, ctx);
    }
    TEMPLATE(T, mat_init) (X, M->c, M->c-rk, ctx);
    if (rk != M->c) 
    {
        /* Invert permutation */
        Qi = flint_malloc(M->c * sizeof(*Qi));
        for (i = 0; i < M->c; ++i) Qi[Q[i]] = i;

        /* Mssign unit vectors to non-pivot columns */
        for (i = M->c-1; i >= rk; --i) TEMPLATE(T, one) (&X->rows[Qi[i]][i-rk], ctx);
        for (i = rk-1; i >= 0; --i) 
        {
            Urow = &U->rows[i];
            Xrow = X->rows[Qi[i]];
            for (j = 1; j < Urow->nnz; ++j) 
            {
                e = &Urow->entries[j];
                /* Do in-place row elimination */
                TEMPLATE(T, neg) (cc, e->val, ctx);
                if (e->ind < rk) _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (Xrow, X->rows[Qi[e->ind]], X->c, cc, ctx);
                else  TEMPLATE(T, sub) (&Xrow[e->ind-rk], &Xrow[e->ind-rk], e->val, ctx);
            }

        }
        flint_free(Qi);
    }
    flint_free(Q);
    TEMPLATE(T, clear) (cc, ctx);
    TEMPLATE(T, sparse_mat_clear) (U, ctx);
    return X->c;
}

#endif
