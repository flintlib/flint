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

static void heap_up(slong *heap, slong *heap_idx, slong *scores, slong pos) {
    const slong c = heap[pos];
    slong nc, npos;
    for(; pos > 0; pos = npos) {
        npos = (pos-1)/2;
        nc = heap[npos];
        if(scores[c] >= scores[nc]) break;

        heap[pos] = nc;
        heap_idx[nc] = pos;
    }
    heap[pos] = c;
    heap_idx[c] = pos;
}

static void heap_down(slong *heap, slong *heap_idx, slong *scores, const slong size, slong pos) {
    const slong c = heap[pos];
    slong nc, npos;
    for(; pos < (size-1)/2; pos = npos) {
        npos = 2*pos+1;
        if(npos+1 < size && scores[heap[npos]] > scores[heap[npos+1]]) ++npos;
        nc = heap[npos];
        if(scores[c] <= scores[nc]) break;

        heap[pos] = nc;
        heap_idx[nc] = pos;
    }
    heap[pos] = c;
    heap_idx[c] = pos;
}

void nmod_sparse_mat_lu(slong *P, slong *Q, 
                        nmod_sparse_mat_t L, nmod_sparse_mat_t U, 
                        const nmod_sparse_mat_t A) {
    slong i, j, r, c;
    nmod_sparse_mat_t Lt;
    nmod_sparse_mat_init(Lt, A->c, A->r, A->mod);
    nmod_sparse_mat_set(U, A);
    nmod_sparse_mat_transpose(Lt, A);

    // Set up permutations
    slong remr = U->r, remc = Lt->r;
    for(r=0; r<U->r; ++r) 
    {
        if(!U->rows[r].nnz) P[r] = --remr; 
        else P[r] = -1;
    }
    for(c=0; c<Lt->r; ++c) 
    {
        if(!Lt->rows[c].nnz) Q[c] = --remc; 
        else Q[c] = -1;
    }

    /* Make heap of nonzero columns by size */
    slong *heap, *heap_idx, *scores, heap_size = A->c;
    heap = flint_malloc(A->c*sizeof(*heap));
    scores = flint_malloc(A->c*sizeof(*scores));
    heap_idx = flint_malloc(A->c*sizeof(*heap_idx));
    for(i=0; i<A->c; ++i) {
        scores[i] = Lt->rows[i].nnz; // TODO: randomized tiebreaker
        heap[i] = i;
        heap_up(heap, heap_idx, scores, i);
    }

    /* Run elimination */
    slong pc, pr;
    slong numr = 0, numc = 0;
    nmod_sparse_vec_struct *pcol, *prow, *row, *col;
    for(heap_size=A->c; heap_size > 0; --heap_size) {
        /* Get lowest weight column (top of heap) */
        pc = heap[0];
        pcol = &Lt->rows[pc];
        heap[0] = heap[heap_size-1];
        heap_down(heap, heap_idx, scores, heap_size, 0);
        if(pcol->nnz==0) continue; // Empty columns already dealt with
        Q[pc] = numc++; // Move pivot column to front
        
        /* Get lowest weight incident row */
        pr = pcol->entries[0].ind, prow = &U->rows[pr];
        for(j=1; j<pcol->nnz; ++j) {
            r = pcol->entries[j].ind, row = &U->rows[r];
            if(row->nnz < prow->nnz) pr = r, prow = row;
        }
        P[pr] = numr++; // Move pivot row to front

        /* Invert pivot */
        mp_limb_t cinv = nmod_inv(nmod_sparse_vec_at(prow, c), A->mod);

        /* Gaussian eliminate rows */
        mp_limb_t cc;
        for(j=0; j<pcol->nnz; ++j) {
            r = pcol->entries[j].ind, row = &U->rows[r];
            if(P[r] >= 0) continue; // Skip previous pivot rows
            cc = nmod_mul(cinv, nmod_sparse_vec_at(row, c), A->mod);
            nmod_sparse_vec_scalar_addmul(row, row, prow, cc, A->mod);
            if(row->nnz==0) P[r] = --remr;
        }

        /* Gaussian eliminate cols */
        for(j=0; j<prow->nnz; ++j) {
            c = prow->entries[j].ind, col = &Lt->rows[c];
            if(Q[c] >= 0) continue; // Skip previous pivot columns
            cc = nmod_mul(cinv, nmod_sparse_vec_at(col, c), A->mod);
            nmod_sparse_vec_scalar_addmul(col, col, pcol, cc, A->mod);
            if(col->nnz==0) Q[c] = --remc;
        }
    }
    /* Reorder cols in U and L^t */
    nmod_sparse_mat_permute_cols(U, Q);
    nmod_sparse_mat_permute_cols(Lt, P);

    /* Transpose L^t */
    nmod_sparse_mat_transpose(L, Lt);
    nmod_sparse_mat_clear(Lt);
}
