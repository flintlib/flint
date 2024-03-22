/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "gr_sparse_mat.h"

int gr_lil_mat_transpose(gr_lil_mat_t B, const gr_lil_mat_t A, gr_ctx_t ctx)
{
    slong r, c, i, j, nz_idx, sz;
    gr_sparse_vec_struct *Arow, *Brow;
    int status = GR_SUCCESS;

    sz = ctx->sizeof_elem;
    r = gr_sparse_mat_nrows(A, ctx);
    c = gr_sparse_mat_ncols(A, ctx);

    if (r != gr_sparse_mat_ncols(B, ctx) || c != gr_sparse_mat_nrows(B, ctx))
        return GR_DOMAIN;

    // TODO: handle
    if (A == B)
        return GR_DOMAIN;
    
    /* Get number of nnzs in each column of A (thus each row of B) */
    for (j = 0; j < c; ++j) 
    {
        B->rows[j].nnz = 0;
    }
    for (i = 0; i < A->r; ++i) 
    {
        Arow = &A->rows[i];
        for (nz_idx = 0; nz_idx < A->rows[i].nnz; ++nz_idx) 
        {
            B->rows[Arow->inds[nz_idx]].nnz += 1;
        }
    }
    /* Allocate space for nnz and reset counters */
    for (j = 0; j < c; ++j) 
    {
        Brow = &B->rows[j];
        gr_sparse_vec_fit_nnz(Brow, Brow->nnz, ctx);
        Brow->nnz = 0;
    }
    /* Put entries into transposed matrix */
    for (i = 0; i < r; ++i)
    {
        Arow = &A->rows[i];
        for (nz_idx = 0; nz_idx < Arow->nnz; ++nz_idx) 
        {
            Brow = &B->rows[Arow->inds[nz_idx]];
            Brow->inds[Brow->nnz] = i;
            status |= gr_set(GR_ENTRY(Brow->nzs, Brow->nnz, sz), GR_ENTRY(Arow->nzs, nz_idx, sz), ctx);
            (Brow->nnz)++;   
        }
    }
    return status;
}
