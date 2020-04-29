/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"
#include "hashmap.h"
#include "heap.h"
#include "perm.h"
#include "longlong.h"

static void update_heap(heap_t h, const fmpz_sparse_mat_with_transpose_t MT, slong r, slong *osupp, slong onnz)
{
    slong i = 0, j = 0, oc, nc, c;
    fmpz_sparse_vec_struct *row = &MT->M->rows[r];
    while (1)
    {
        oc = (i==onnz) ? MT->M->c : osupp[i];
        nc = (j==row->nnz) ? MT->M->c : row->entries[j].ind;
        c = FLINT_MIN(oc, nc);
        if (c >= MT->M->c) break;
        if (oc != nc) heap_adjust(h, c, MT->cols[c].num);
        if (c==oc) ++i;
        if (c==nc) ++j;
    }
}

slong
fmpz_sparse_mat_fflu(fmpz *D, slong *P, slong *Q, fmpz_sparse_mat_t L, fmpz_sparse_mat_t U, 
                            const fmpz_sparse_mat_t M)
{
    slong i, j, r, c, pr, pc, rank, remr, remc, *supp, nnz;
    heap_t h;
    fmpz_t pivot, one, tmp;
    fmpz_sparse_vec_struct *prow, *row;
    hashmap_struct *hcol;
    fmpz_t *hval;
    fmpz_sparse_mat_with_transpose_t UT;

    for (i = 0; i < M->r; ++i) fmpz_one(&D[i]);
    if (M->r == 0 || M->c == 0) 
    {
        fmpz_sparse_mat_zero (L);
        fmpz_sparse_mat_zero (U);
        for (i = 0; i < M->r; ++i) P[i] = i;
        for (i = 0; i < M->c; ++i) Q[i] = i;
        return 0;
    }
    fmpz_init(pivot);
    fmpz_init(tmp);
    fmpz_init_set_ui(one, UWORD(1));
    fmpz_sparse_mat_zero(L);
    fmpz_sparse_mat_set (U, M);
    _fmpz_sparse_mat_with_transpose_init(UT, U);

    /* Initialize permutations */
    remr = M->r, remc = M->c;
    for (r = 0; r < M->r; ++r) 
    {
        if (!U->rows[r].nnz) P[r] = --remr; 
        else P[r] = -1;
    }
    for (c = 0; c < M->c; ++c) 
    {
        if (!UT->cols[c].num) Q[c] = --remc; 
        else Q[c] = -1;
    }
    
    /* Make heap of nonzero columns by size */
    heap_init(h, M->c);
    for (c = 0; c < M->c; ++c)
        heap_push(h, UT->cols[c].num);

    /* Run elimination */
    rank = 0;

    while (h->num > 0)
    {
        /* Get lowest weight column (top of heap) */
        pc = heap_pop(h, NULL);
        hcol = &UT->cols[pc];
        
        /* Get lowest weight incident row */
        prow = NULL, pr = -1;
        for (i = 0; i < hcol->num; ++i)
        {
            r = hcol->keys[i];
            row = &U->rows[r];
            if (prow == NULL || row->nnz < prow->nnz)
                pr = r, prow = row;
        }
        if (pr == -1) {if (Q[pc] == -1) Q[pc] = --remc; continue;}
        P[pr] = Q[pc] = rank++; /* Move pivot row and col to front */

        /* Get pivot */
        fmpz_set(pivot, *fmpz_sparse_vec_at(prow, pc));
        fmpz_sparse_vec_set_entry(&L->rows[pr], pc, one);

        /* Remove pivot row from transpose */
        for (j = 0; j < prow->nnz; ++j)
            hashmap_rem(&UT->cols[prow->entries[j].ind], pr);

        /* Use pivot row to fraction-free Gaussian eliminate other incident rows */
        while (hcol->num > 0)
        {
            r = hcol->keys[0];
            nnz = _fmpz_sparse_vec_support(&supp, &U->rows[r]);
            fmpz_mul(&D[r], &D[r], pivot);
            hval = hcol->vals[0], fmpz_set(tmp, *hval);
            fmpz_sparse_vec_scalar_mul_fmpz(&U->rows[r], &U->rows[r], pivot);
            fmpz_sparse_vec_scalar_mul_fmpz(&L->rows[r], &L->rows[r], pivot);
            fmpz_sparse_vec_set_entry(&L->rows[r], pc, tmp);
            _fmpz_sparse_mat_with_transpose_gauss_elim_col(UT, pr, r, pc);
            update_heap(h, UT, r, supp, nnz);
            if (U->rows[r].nnz == 0) P[r] = --remr;
        }
    }
    heap_clear(h);
    _fmpz_sparse_mat_with_transpose_clear(UT);

    /* Reorder rows and cols in L and U */
    fmpz_sparse_mat_permute_cols (L, Q);
    fmpz_sparse_mat_permute_rows (L, P);
    fmpz_sparse_mat_permute_cols (U, Q);
    fmpz_sparse_mat_permute_rows (U, P);
    fmpz_clear(pivot);
    fmpz_clear(tmp);
    return rank;
}
