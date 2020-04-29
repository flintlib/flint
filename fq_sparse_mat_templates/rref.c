/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include <string.h>
#include "templates.h"

slong TEMPLATE(T, sparse_mat_rref) (TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx)
{
    slong *P;
    slong j, r, c, pr, pc, rank, remr;
    TEMPLATE(T, t) cinv, cc;
    TEMPLATE(T, sparse_mat_t) Mt;
    TEMPLATE(T, sparse_vec_struct) *pcol, *prow, *row, *col;

    if (M->r == 0 || M->c == 0) return 0;
    P = flint_malloc(M->r*sizeof(*P));
    TEMPLATE(T, init) (cinv, ctx);
    TEMPLATE(T, init) (cc, ctx);
    TEMPLATE(T, sparse_mat_init) (Mt, M->c, M->r, ctx);

    /* Make transpose for doing column eliminations */
    TEMPLATE(T, sparse_mat_transpose) (Mt, M, ctx);
    
    /* Set up permutations */
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

        TEMPLATE(T, inv) (cinv, *TEMPLATE(T, sparse_vec_at) (prow, pc, ctx), ctx);
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (prow, prow, cinv, ctx);

        /* Gaussian eliminate rows */
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->entries[j].ind, row = &M->rows[r];
            if (r == pr) {TEMPLATE(T, zero) (pcol->entries[j].val, ctx); continue;}

            TEMPLATE(T, neg) (cc, *TEMPLATE(T, sparse_vec_at) (row, pc, ctx), ctx);
            TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (row, row, prow, cc, ctx);
            if (row->nnz == 0 || row->entries[0].ind >= M->c) P[r] = --remr;
        }
        /* Gaussian eliminate cols */
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (pcol, pcol, cinv, ctx);
        for (j = 0; j < prow->nnz; ++j)
        {
            c = prow->entries[j].ind, col = &Mt->rows[c];
            if (c >= M->c || c == pc) continue;
            TEMPLATE(T, neg) (cc, *TEMPLATE(T, sparse_vec_at) (col, pr, ctx), ctx);
            TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (col, col, pcol, cc, ctx);
        }
        rank += 1;
    }
    /* Reorder rows */
    TEMPLATE(T, sparse_mat_permute_rows) (M, P, ctx);

    flint_free(P);
    TEMPLATE(T, clear) (cinv, ctx);
    TEMPLATE(T, clear) (cc, ctx);
    TEMPLATE(T, sparse_mat_clear) (Mt, ctx);

    return rank;
}

#endif
