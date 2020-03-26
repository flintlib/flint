/*
    Copyright (C) 2011 Fredrik Johansson

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

void
nmod_sparse_mat_transpose(nmod_sparse_mat_t B, const nmod_sparse_mat_t A)
{

    slong i, j;
    B->c = A->r;
    B->nnz = A->nnz;
    B->entries = flint_realloc(B->entries, B->nnz * sizeof(*B->entries));
    memset(B->row_starts, 0, B->r * sizeof(*B->row_starts));
    memset(B->row_nnz, 0, B->r * sizeof(*B->row_nnz));
    if(A->nnz == 0) return;

    /* Clear row counts for B and set row starts */
    nmod_sparse_mat_entry_struct *Ae, *Be;
    Ae = A->entries;
    for(i=0; i<A->nnz; ++i)
        if(Ae[i].col + 1 < B->r)
            B->row_starts[Ae[i].col + 1] += 1;
    for(i=1; i<B->r; ++i)
        B->row_starts[i] += B->row_starts[i-1];

    /* Assign entries */
    Ae = A->entries;
    Be = B->entries;
    for(i=0; i<A->r; ++i) {
        for(j=0; j<A->row_nnz[i]; ++j, ++Ae) {
            slong pos = B->row_starts[Ae->col] + B->row_nnz[Ae->col]++;
            Be[pos].col = i;
            Be[pos].val = Ae->val;
        }
        flint_printf("\n");
    }
}
