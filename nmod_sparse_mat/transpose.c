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
    slong r, c, i, j;
    nmod_sparse_entry_struct *Ae, *Be;
    nmod_sparse_vec_struct *row;

    /* Get number of nnzs in each column of A (thus each row of B) */
    for (c = 0; c < A->c; ++c) 
    {
        B->rows[c].nnz = 0;
    }
    for (r = 0; r < A->r; ++r) 
    {
        for (i = 0; i < A->rows[r].nnz; ++i) 
        {
            c = A->rows[r].entries[i].ind;
            if (c >= A->c) break;
            B->rows[c].nnz += 1;
        }
    }
    /* Allocate space for nnz and reset counters */
    for (c = 0; c < A->c; ++c) 
    {
        row = &B->rows[c];
        if (row->nnz == 0) nmod_sparse_vec_clear(row);
        else row->entries = flint_realloc(row->entries, row->nnz*sizeof(*row->entries));
        row->nnz = 0;
    }
    /* Put entries into transposed matrix */
    for (r = 0; r < A->r; ++r)
    {
        for (i = 0; i < A->rows[r].nnz; ++i) 
        {
            Ae = &A->rows[r].entries[i];
            c = Ae->ind;
            if (c >= A->c) break;
            j = B->rows[c].nnz++;
            Be = &B->rows[c].entries[j];
            Be->ind = r, Be->val = Ae->val;
        }
    }
    B->c_off = 0;
}
