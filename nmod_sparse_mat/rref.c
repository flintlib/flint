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

slong nmod_sparse_mat_rref(nmod_sparse_mat_t A)
{
    if (A->r == 0 || A->c == 0) return 0;
    slong *P;
    slong i, j, r, c, pr, pc, rank, remr;
    nmod_sparse_mat_t At;
    nmod_sparse_vec_struct *pcol, *prow, *row, *col;
    mp_limb_t cinv, cc;

    nmod_sparse_mat_init(At, A->c, A->r, A->mod);
    nmod_sparse_mat_transpose(At, A);
    
    /* Set up permutations */
    P = flint_malloc(A->r*sizeof(*P));
    remr = A->r;
    for (r = 0; r<A->r; ++r) 
    {
        if (!A->rows[r].nnz || A->rows[r].entries[0].ind >= A->c) P[r] = --remr; 
        else P[r] = -1;
    }
    
    /* Run elimination */
    rank = 0;
    for (pc=0; pc<A->c; ++pc)
    {
        pcol = &At->rows[pc];

        /* Get lowest weight incident row not used as previous pivot */
        pr = -1, prow = NULL;
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->entries[j].ind, row = &A->rows[r];
            if (P[r] >= 0) continue;
            if (pr==-1 || (row->nnz < prow->nnz)) pr = r, prow = row;
        }
        if(pr == -1) continue;
        P[pr] = rank; 

        cinv = nmod_inv(nmod_sparse_vec_at(prow, pc), A->mod);
        nmod_sparse_vec_scalar_mul(prow, prow, cinv, A->mod);

        /* Gaussian eliminate rows */
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->entries[j].ind, row = &A->rows[r];
            if(r==pr) {pcol->entries[j].val = UWORD(0); continue;}

            cc = nmod_neg(nmod_sparse_vec_at(row, pc), A->mod);
            nmod_sparse_vec_scalar_addmul(row, row, prow, cc, A->mod);
            if (row->nnz == 0 || row->entries[0].ind >= A->c) P[r] = --remr;
        }
        /* Gaussian eliminate cols */
        nmod_sparse_vec_scalar_mul(pcol, pcol, cinv, A->mod);
        for (j = 0; j < prow->nnz; ++j)
        {
            c = prow->entries[j].ind, col = &At->rows[c];
            if(c >= A->c || c==pc) continue;
            cc = nmod_neg(nmod_sparse_vec_at(col, pr), A->mod);
            nmod_sparse_vec_scalar_addmul(col, col, pcol, cc, A->mod);
        }
        rank += 1;
    }
    nmod_sparse_mat_clear(At);

    /* Reorder rows */
    nmod_sparse_mat_permute_rows(A, P);
    flint_free(P);
    return rank;
}
