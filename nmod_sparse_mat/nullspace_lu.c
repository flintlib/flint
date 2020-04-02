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

slong nmod_sparse_mat_nullspace_lu(nmod_mat_t X, const nmod_sparse_mat_t A)
{
    int good = 1;
    slong rk, *P, *Q, *Qi, i, j;
    nmod_sparse_mat_t L, U;
    nmod_sparse_entry_struct *e;
    nmod_sparse_vec_struct *Urow;
    mp_limb_t *Xrow;
    P = flint_malloc(A->r * sizeof(*P));
    Q = flint_malloc(A->c * sizeof(*Q));
    nmod_sparse_mat_init(L, A->r, A->c, A->mod);
    nmod_sparse_mat_init(U, A->r, A->c, A->mod);
    rk = nmod_sparse_mat_lu(P, Q, L, U, A);
    flint_free(P);
    nmod_sparse_mat_clear(L);
    for(i=0; i<rk; ++i)
        nmod_sparse_vec_scalar_mul(&U->rows[i], &U->rows[i], nmod_inv(U->rows[i].entries[0].val, A->mod), A->mod);
    nmod_mat_init(X, A->c, A->c-rk, A->mod.n);
    if (rk != A->c) 
    {
        /* Invert permutation */
        Qi = flint_malloc(A->c * sizeof(*Qi));
        for (i = 0; i < A->c; ++i) Qi[Q[i]] = i;

        /* Assign unit vectors to non-pivot columns */
        for (i = A->c-1; i >= rk; --i) X->rows[Qi[i]][i-rk] = 1;
        for (i = rk-1; i >= 0; --i) {
            Urow = &U->rows[i];
            Xrow = X->rows[Qi[i]];
            for(j = 1; j<Urow->nnz; ++j) {
                e = &Urow->entries[j];
                /* Do in-place row elimination */
                if(e->ind < rk) _nmod_vec_scalar_addmul_nmod(Xrow, X->rows[Qi[e->ind]], X->c, nmod_neg(e->val, A->mod), A->mod);
                else Xrow[e->ind-rk] = nmod_sub(Xrow[e->ind-rk], e->val, A->mod);
            }

        }
        flint_free(Qi);
    }
    flint_free(Q);
    nmod_sparse_mat_clear(U);
    return X->c;
}
