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
#include "nmod_vec.h"

void
nmod_sparse_mat_add(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B)
{
    slong i;

    C->entries = flint_realloc(C->entries, (A->nnz + B->nnz)*sizeof(*C->entries));
    memset(C->entries, 0, (A->nnz + B->nnz)*sizeof(*C->entries));
    C->nnz = 0;

    for (i = 0; i < C->r; i++)
    {
        C->row_starts[i] = C->nnz;
        nmod_sparse_mat_entry_struct *Ae = A->entries + A->row_starts[i];
        nmod_sparse_mat_entry_struct *Be = B->entries + B->row_starts[i];
        nmod_sparse_mat_entry_struct *Ce = C->entries + C->row_starts[i];
        slong j = 0;
        slong k = 0;
        slong nnz = 0;
        /* Interleave ith rows until one runs out */
        while(j < A->row_nnz[i] && k < B->row_nnz[i]) {
            slong col = Ce[nnz].col = FLINT_MIN(Ae[j].col, Be[k].col);
            if(Ae[j].col == col) Ce[nnz].val = Ae[j++].val;
            if(Be[k].col == col) Ce[nnz].val = nmod_add(Ce[nnz].val, Be[k++].val, C->mod);
            if(Ce[nnz].val != UWORD(0)) ++nnz;
        }
        /* Add remainder of A row */
        for(; j<A->row_nnz[i]; ++j, ++nnz) Ce[nnz] = Ae[j];
        /* Add remainder of B row */
        for(; k<B->row_nnz[i]; ++k, ++nnz) Ce[nnz]= Be[k];

        C->row_nnz[i] = nnz;
        C->nnz += nnz;
    }
    _nmod_sparse_mat_set_c(C);
    C->entries = realloc(C->entries, C->nnz*sizeof(*C->entries));
}
