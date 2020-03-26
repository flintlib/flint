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
void nmod_sparse_mat_mul_vec(mp_ptr y, const nmod_sparse_mat_t A, const mp_ptr x) {
    slong i,j;
    memset(y, 0, A->r * sizeof(*y));
    nmod_sparse_mat_entry_struct *e = A->entries;
    for(i=0; i<A->r; ++i) {
        for(j=0; j<A->row_nnz[i]; ++j, ++e) {
            if(e->val==1) 
                y[i] = nmod_add(y[i], x[e->col], A->mod);
            else if(e->val==A->mod.n-1) 
                y[i] = nmod_sub(y[i], x[e->col], A->mod);
            else 
                y[i] = nmod_add(y[i], nmod_mul(x[e->col], e->val, A->mod), A->mod);
        }
    }
}

