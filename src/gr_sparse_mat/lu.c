/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "gr_sparse_mat.h"

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
int gr_lil_mat_lu(
    slong *res_rank, slong *P, slong *Q, 
    gr_lil_mat_t L, gr_lil_mat_t U, 
    const gr_lil_mat_t A, gr_ctx_t ctx
)
{
    slong i, j, r, c, rank, pr, pc, remr, remc;
    slong *heap, *heap_idx, *scores, heap_size;
    gr_lil_mat_t Lt;
    gr_sparse_vec_struct *pcol, *prow, *row, *col;
    gr_ptr cinv, cc;
    int status = GR_SUCCESS;

    if (A->r == 0 || A->c == 0 || A->nnz == 0) 
    {
        *res_rank = 0;
        gr_lil_mat_zero(L, ctx);
        gr_lil_mat_zero(U, ctx);
        for (i = 0; i < A->r; ++i) P[i] = i;
        for (j = 0; j < A->c; ++j) Q[j] = j;
        return GR_SUCCESS;
    }
    GR_TMP_INIT2(cinv, cc, ctx);
    gr_lil_mat_init(Lt, L->c, L->r, ctx);
    status |= gr_lil_mat_transpose(Lt, A, ctx);
    status |= gr_lil_mat_set(U, A, ctx);
    
    /* Set up permutations */
    remr = A->r, remc = A->c;
    for (r = 0; r < A->r; ++r) 
    {
        if (!U->rows[r].nnz) P[r] = --remr; 
        else P[r] = -1;
    }
    if (Q != NULL)
    {
        for (j = 0; j < A->c; ++j) 
        {
            if (!Lt->rows[j].nnz) Q[j] = --remc; 
            else Q[j] = -1;
        }
        /* Make heap of nonzero columns by size */
        heap_size = A->c;
        heap = flint_malloc(A->c*sizeof(*heap));
        scores = flint_malloc(A->c*sizeof(*scores));
        heap_idx = flint_malloc(A->c*sizeof(*heap_idx));
        for (j = 0; j < A->c; ++j)
        {
            scores[j] = Lt->rows[j].nnz; /* TODO: randomized tiebreaker */
            heap[j] = j;
            heap_up(heap, heap_idx, scores, j);
        }
    }
    
    /* Run elimination */
    rank = 0;
    for (heap_size = A->c; heap_size > 0; )
    {

        /* Get pivot column */
        if (Q != NULL)
        {
            /* Get lowest weight column (top of heap) */
            pc = heap[0];
            pcol = &Lt->rows[pc];
            heap[0] = heap[--heap_size];
            heap_down(heap, heap_idx, scores, heap_size, 0);
            if (pcol->nnz == 0) continue; /* Empty columns already dealt with */
            Q[pc] = rank; /* Move pivot column to front */
        }
        else
        {
            pc = A->c - heap_size--;
            pcol = &Lt->rows[pc];
            if (pcol->nnz == 0) continue; /* Nothing to do */
        }
        
        /* Get lowest weight incident row */
        pr = pcol->inds[0], prow = &U->rows[pr];
        for (j = 1; j < pcol->nnz; ++j)
        {
            r = pcol->inds[j], row = &U->rows[r];
            if (row->nnz < prow->nnz) pr = r, prow = row;
        }
        P[pr] = rank; /* Move pivot row to front */
        
        /* Invert pivot */
        status |= gr_inv(cinv, gr_sparse_vec_find_entry(prow, pc, ctx), ctx);

        /* Lower triangular matrix will have ones on the diagonal */
        status |= gr_sparse_vec_mul_scalar(pcol, pcol, cinv, ctx);

        /* Gaussian eliminate lower rows in U incident on pivot column */
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->inds[j], row = &U->rows[r];
            if (P[r] >= 0) continue; /* Skip previous pivot rows */
            
            status |= gr_mul(cc, cinv, gr_sparse_vec_find_entry(row, pc, ctx), ctx);
            status |= gr_neg(cc, cc, ctx);
            status |= gr_sparse_vec_addmul_scalar(row, prow, cc, ctx);
            if (row->nnz == 0) P[r] = --remr;
        }
        /* Gaussian eliminate lower cols in L incident on pivot row */
        for (j = 0; j < prow->nnz; ++j)
        {
            c = prow->inds[j], col = &Lt->rows[c];
            if ((Q == NULL && (c >= A->c || c<=pc)) || (Q != NULL && Q[c] >= 0))
                continue; /* Skip previous pivot columns */
            status |= gr_neg(cc, gr_sparse_vec_find_entry(col, pr, ctx), ctx);
            status |= gr_sparse_vec_addmul_scalar(col, pcol, cc, ctx);
            
            if (Q != NULL)
            {
                if (col->nnz == 0) Q[c] = --remc;
                scores[c] = col->nnz;
                heap_up(heap, heap_idx, scores, heap_idx[c]);
                heap_down(heap, heap_idx, scores, heap_size, heap_idx[c]);
            }
        }
        rank += 1;
    }
    /* Fix nnz */
    Lt->nnz = 0;
    for (j = 0; j < A->c; ++j)
    {
        Lt->nnz += Lt->rows[j].nnz;
    }
    U->nnz = 0;
    for (j = 0; j < A->r; ++j)
    {
        U->nnz += U->rows[j].nnz;
    }

    /* Transpose L^t */
    status |= gr_lil_mat_transpose(L, Lt, ctx);

    /* Reorder rows and cols in L and U */
    status |= gr_lil_mat_permute_rows(L, P, ctx);
    status |= gr_lil_mat_permute_rows(U, P, ctx);
    if (Q != NULL)
    {
        status |= gr_lil_mat_permute_cols(L, Q, ctx);
        status |= gr_lil_mat_permute_cols(U, Q, ctx);
    }
    *res_rank = rank;

    gr_lil_mat_clear(Lt, ctx);
    GR_TMP_INIT2(cinv, cc, ctx);
    return status;
}
