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

slong nmod_sparse_mat_rref(nmod_sparse_mat_t M)
{
    slong *P;
    slong j, r, c, pr, pc, rank, remr;
    nmod_sparse_mat_t Mt;
    nmod_sparse_vec_struct *pcol, *prow, *row, *col;
    mp_limb_t cinv, cc;

    if (M->r == 0 || M->c == 0) return 0;
    nmod_sparse_mat_init(Mt, M->c, M->r, M->mod);
    nmod_sparse_mat_transpose(Mt, M);
    
    /* Set up permutations */
    P = flint_malloc(M->r*sizeof(*P));
    remr = M->r;
    for (r = 0; r < M->r; ++r) 
    {
        if (!M->rows[r].nnz || M->rows[r].entries[0].ind >= M->c) P[r] = --remr; 
        else P[r] = -1;
    }
    
    /* Run elimination */
    rank = 0;
    for (pc = 0; pc < M->c; ++pc)
    {
        pcol = &Mt->rows[pc];

        /* Get lowest weight incident row not used as previous pivot */
        pr = -1, prow = NULL;
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->entries[j].ind, row = &M->rows[r];
            if (P[r] >= 0) continue;
            if (pr==-1 || (row->nnz < prow->nnz)) pr = r, prow = row;
        }
        if (pr == -1) continue;
        P[pr] = rank; 

        cinv = nmod_inv(*nmod_sparse_vec_at(prow, pc), M->mod);
        nmod_sparse_vec_scalar_mul_nmod(prow, prow, cinv, M->mod);

        /* Gaussian eliminate rows */
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->entries[j].ind, row = &M->rows[r];
            if (r==pr) {pcol->entries[j].val = UWORD(0); continue;}

            cc = nmod_neg(*nmod_sparse_vec_at(row, pc), M->mod);
            nmod_sparse_vec_scalar_addmul_nmod(row, row, prow, cc, M->mod);
            if (row->nnz == 0 || row->entries[0].ind >= M->c) P[r] = --remr;
        }
        /* Gaussian eliminate cols */
        nmod_sparse_vec_scalar_mul_nmod(pcol, pcol, cinv, M->mod);
        for (j = 0; j < prow->nnz; ++j)
        {
            c = prow->entries[j].ind, col = &Mt->rows[c];
            if (c >= M->c || c==pc) continue;
            cc = nmod_neg(*nmod_sparse_vec_at(col, pr), M->mod);
            nmod_sparse_vec_scalar_addmul_nmod(col, col, pcol, cc, M->mod);
        }
        rank += 1;
    }
    nmod_sparse_mat_clear(Mt);

    /* Reorder rows */
    nmod_sparse_mat_permute_rows(M, P);
    flint_free(P);
    return rank;
}
