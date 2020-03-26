/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

/* TODO: save reductions to end? */
void nmod_sparse_mat_mul_mat(nmod_mat_t Y, const nmod_sparse_mat_t A, const nmod_mat_t X) {
    slong i,j,k;
    if(Y->r == 0 || Y->c == 0) return;
    memset(Y->entries, 0, Y->r * Y->c * sizeof(*Y->entries));
    nmod_sparse_mat_entry_struct *e = A->entries;
    for(i=0; i<A->r; ++i) {
        for(j=0; j<A->row_nnz[i]; ++j, ++e) {
            if(e->val==1) 
                _nmod_vec_add(Y->rows[i], Y->rows[i], X->rows[e->col], X->c, A->mod);
            else if(e->val==A->mod.n-1) 
                _nmod_vec_sub(Y->rows[i], Y->rows[i], X->rows[e->col], X->c, A->mod);
            else 
                _nmod_vec_scalar_addmul_nmod(Y->rows[i], X->rows[e->col], X->c, e->val, A->mod);
        }
    }
}
