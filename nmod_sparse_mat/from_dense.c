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

/* Assumes B->c == A->r */
void
nmod_sparse_mat_from_dense(nmod_sparse_mat_t B, const nmod_mat_t A)
{
    slong i, j;
    if(A->r == 0) return;
    if(A->c == 0) {nmod_sparse_mat_zero(B); return;}
    B->entries = flint_realloc(B->entries, A->r * A->c * sizeof(*B->entries));
    B->nnz = 0;
    B->c = 0;
    for(i=0; i < A->r; ++i) {
        B->row_starts[i] = B->nnz;
        B->row_nnz[i] = 0;
        for(j=0; j < A->c; ++j) {
            if(A->rows[i][j] != UWORD(0)) {
                B->entries[B->nnz].col = j;
                B->entries[B->nnz].val = A->rows[i][j];
                B->row_nnz[i]++;
                B->nnz++;
            }
        }
    }
    _nmod_sparse_mat_set_c(B);
    if(B->nnz > 0)
        B->entries = flint_realloc(B->entries, B->nnz * sizeof(*B->entries));
    else {
        flint_free(B->entries);
        B->entries = NULL;
    }
}
