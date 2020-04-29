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

slong nmod_sparse_mat_nullspace_lu(nmod_mat_t X, const nmod_sparse_mat_t M)
{
    slong rk, *P, *Q, *Qi, i, j;
    nmod_sparse_mat_t L, U;
    nmod_sparse_entry_struct *e;
    nmod_sparse_vec_struct *Urow;
    mp_limb_t *Xrow;
    P = flint_malloc(M->r * sizeof(*P));
    Q = flint_malloc(M->c * sizeof(*Q));
    nmod_sparse_mat_init(L, M->r, M->c, M->mod);
    nmod_sparse_mat_init(U, M->r, M->c, M->mod);
    rk = nmod_sparse_mat_lu(P, Q, L, U, M);
    flint_free(P);
    nmod_sparse_mat_clear(L);
    for (i = 0; i < rk; ++i)
        nmod_sparse_vec_scalar_mul_nmod(&U->rows[i], &U->rows[i], nmod_inv(U->rows[i].entries[0].val, M->mod), M->mod);
    nmod_mat_init(X, M->c, M->c-rk, M->mod.n);
    if (rk != M->c) 
    {
        /* Invert permutation */
        Qi = flint_malloc(M->c * sizeof(*Qi));
        for (i = 0; i < M->c; ++i) Qi[Q[i]] = i;

        /* Mssign unit vectors to non-pivot columns */
        for (i = M->c-1; i >= rk; --i) X->rows[Qi[i]][i-rk] = 1;
        for (i = rk-1; i >= 0; --i) 
        {
            Urow = &U->rows[i];
            Xrow = X->rows[Qi[i]];
            for (j = 1; j < Urow->nnz; ++j) 
            {
                e = &Urow->entries[j];
                /* Do in-place row elimination */
                if (e->ind < rk) _nmod_vec_scalar_addmul_nmod(Xrow, X->rows[Qi[e->ind]], X->c, nmod_neg(e->val, M->mod), M->mod);
                else Xrow[e->ind-rk] = nmod_sub(Xrow[e->ind-rk], e->val, M->mod);
            }

        }
        flint_free(Qi);
    }
    flint_free(Q);
    nmod_sparse_mat_clear(U);
    return X->c;
}
