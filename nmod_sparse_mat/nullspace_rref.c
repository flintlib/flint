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

slong nmod_sparse_mat_nullspace_rref(nmod_mat_t X, const nmod_sparse_mat_t A)
{
    slong i, j, rk, numc, *Q;
    nmod_sparse_mat_t R;
    nmod_sparse_vec_struct *Rrow;
    mp_limb_t *Xrow;
    nmod_sparse_entry_struct *le, *e;
    nmod_sparse_mat_init(R, A->r, A->c, A->mod);
    nmod_sparse_mat_set(R, A);
    rk = nmod_sparse_mat_rref(R);
    nmod_mat_init(X, A->c, A->c-rk, A->mod.n);
    if(rk != A->c) 
    {
        numc = 0;
        /* Mark which cols are pivots and enumerate the nonpivots */
        Q = flint_calloc(A->c, sizeof(*Q));
        for (i = 0; i < rk; ++i) 
            Q[R->rows[i].entries->ind] = -1;
        for (i = 0; i < A->c; ++i)
            if(Q[i]==UWORD(0)) Q[i] = numc++, X->rows[i][Q[i]] = 1;

        /* For each pivot col, set the corresponding row in X as */
        /* the negative of the associated row in R (reordered by Q) */
        for (i = 0; i < rk; ++i) 
        {
            Rrow = &R->rows[i];
            Xrow = X->rows[Rrow->entries[0].ind];
            for (j = 1; j < Rrow->nnz; ++j) 
                Xrow[Q[Rrow->entries[j].ind]] = nmod_neg(Rrow->entries[j].val, A->mod);
        }
        flint_free(Q);
    }
    nmod_sparse_mat_clear(R);
    return X->c;
}
