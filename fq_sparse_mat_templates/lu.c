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

static void heap_up(slong *heap, slong *heap_idx, slong *scores, slong pos)
{
    const slong c = heap[pos];
    slong nc, npos;
    for (; pos > 0; pos = npos)
    {
        npos = (pos-1)/2;
        nc = heap[npos];
        if (scores[c] >= scores[nc]) break;

        heap[pos] = nc;
        heap_idx[nc] = pos;
    }
    heap[pos] = c;
    heap_idx[c] = pos;
}

static void heap_down(slong *heap, slong *heap_idx, slong *scores, slong size, slong pos)
{
    const slong c = heap[pos];
    slong nc, npos;
    for (; pos < (size-1)/2; pos = npos)
    {
        npos = 2*pos+1;
        if (npos+1 < size && scores[heap[npos]] > scores[heap[npos+1]]) ++npos;
        nc = heap[npos];
        if (scores[c] <= scores[nc]) break;

        heap[pos] = nc;
        heap_idx[nc] = pos;
    }
    heap[pos] = c;
    heap_idx[c] = pos;
}

/* static void print_heap(slong *heap, slong *scores, slong size) 
{
    slong level, i;
    for (level = 1; level <= size; level<<=1) 
    {
        for (i = level; i <= size && i < 2*level; ++i) 
        {
            flint_printf("%wd:%wd,%wd\t", i-1, heap[i-1], scores[heap[i-1]]);
        }
        flint_printf("\n");
    }
}
 */
slong TEMPLATE(T, sparse_mat_lu)(slong *P, slong *Q, 
                        TEMPLATE(T, sparse_mat_t) L, TEMPLATE(T, sparse_mat_t) U, 
                        const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, r, c, pr, pc, rank, remr, remc;
    slong *heap, *heap_idx, *scores, heap_size;
    TEMPLATE(T, t) cinv, cc;
    TEMPLATE(T, sparse_mat_t) Lt;
    TEMPLATE(T, sparse_vec_struct) *pcol, *prow, *row, *col;

    if (M->r == 0 || M->c == 0) 
    {
        TEMPLATE(T, sparse_mat_zero) (L, ctx);
        TEMPLATE(T, sparse_mat_zero) (U, ctx);
        for (i = 0; i < M->r; ++i) P[i] = i;
        for (i = 0; i < M->c; ++i) Q[i] = i;
        return 0;
    }
    TEMPLATE(T, init) (cinv, ctx);
    TEMPLATE(T, init) (cc, ctx);
    TEMPLATE(T, sparse_mat_init) (Lt, L->c, L->r, ctx);
    TEMPLATE(T, sparse_mat_transpose) (Lt, M, ctx);
    TEMPLATE(T, sparse_mat_set) (U, M, ctx);
    
    /* Set up permutations */
    remr = M->r, remc = M->c;
    for (r = 0; r < M->r; ++r) 
    {
        if (!U->rows[r].nnz) P[r] = --remr; 
        else P[r] = -1;
    }
    for (c = 0; c < M->c; ++c) 
    {
        if (!Lt->rows[c].nnz) Q[c] = --remc; 
        else Q[c] = -1;
    }
    
    /* Make heap of nonzero columns by size */
    heap_size = M->c;
    heap = flint_malloc(M->c*sizeof(*heap));
    scores = flint_malloc(M->c*sizeof(*scores));
    heap_idx = flint_malloc(M->c*sizeof(*heap_idx));
    for (i = 0; i < M->c; ++i)
    {
        scores[i] = Lt->rows[i].nnz; /* TODO: randomized tiebreaker */
        heap[i] = i;
        heap_up(heap, heap_idx, scores, i);
    }
    /* Run elimination */
    rank = 0;
    for (heap_size = M->c; heap_size > 0; )
    {
        /* Get lowest weight column (top of heap) */
        pc = heap[0];
        pcol = &Lt->rows[pc];
        heap[0] = heap[--heap_size];
        heap_down(heap, heap_idx, scores, heap_size, 0);
        if (pcol->nnz == 0) continue; /* Empty columns already dealt with */
        Q[pc] = rank; /* Move pivot column to front */
        
        /* Get lowest weight incident row */
        pr = pcol->entries[0].ind, prow = &U->rows[pr];
        for (j = 1; j < pcol->nnz; ++j)
        {
            r = pcol->entries[j].ind, row = &U->rows[r];
            if (row->nnz < prow->nnz) pr = r, prow = row;
        }
        P[pr] = rank; /* Move pivot row to front */
        
        /* Invert pivot */
        TEMPLATE(T, inv) (cinv, *TEMPLATE(T, sparse_vec_at) (prow, pc, ctx), ctx);

        /* Gaussian eliminate rows */
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->entries[j].ind, row = &U->rows[r];
            if (P[r] >= 0) continue; /* Skip previous pivot rows */
            TEMPLATE(T, mul) (cc, cinv, *TEMPLATE(T, sparse_vec_at) (row, pc, ctx), ctx);
            TEMPLATE(T, neg) (cc, cc, ctx);
            TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (row, row, prow, cc, ctx);
            if (row->nnz == 0) P[r] = --remr;
        }
        /* Gaussian eliminate cols */
        TEMPLATE(T, TEMPLATE(sparse_vec_scalar_mul, T)) (pcol, pcol, cinv, ctx);
        for (j = 0; j < prow->nnz; ++j)
        {
            c = prow->entries[j].ind, col = &Lt->rows[c];
            if (Q[c] >= 0) continue; /* Skip previous pivot columns */
            TEMPLATE(T, neg) (cc, *TEMPLATE(T, sparse_vec_at) (col, pr, ctx), ctx);
            TEMPLATE(T, TEMPLATE(sparse_vec_scalar_addmul, T)) (col, col, pcol, cc, ctx);
            if (col->nnz == 0) Q[c] = --remc;
            scores[c] = col->nnz;
            heap_up(heap, heap_idx, scores, heap_idx[c]);
            heap_down(heap, heap_idx, scores, heap_size, heap_idx[c]);

        }
        rank += 1;
    }
    flint_free(heap);
    flint_free(scores);
    flint_free(heap_idx);
    
    TEMPLATE(T, clear) (cinv, ctx);
    TEMPLATE(T, clear) (cc, ctx);

    /* Transpose L^t */
    TEMPLATE(T, sparse_mat_transpose) (L, Lt, ctx);
    TEMPLATE(T, sparse_mat_clear) (Lt, ctx);

    /* Reorder rows and cols in L and U */
    TEMPLATE(T, sparse_mat_permute_cols) (L, Q, ctx);
    TEMPLATE(T, sparse_mat_permute_rows) (L, P, ctx);
    TEMPLATE(T, sparse_mat_permute_cols) (U, Q, ctx);
    TEMPLATE(T, sparse_mat_permute_rows) (U, P, ctx);
    return rank;
}

#endif
