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
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "gr_sparse_vec.h"
#include "gr_sparse_mat.h"

int gr_lil_mat_rref(slong *res_rank, gr_lil_mat_t R, gr_lil_mat_t A, gr_ctx_t ctx)
{
    slong *P;
    slong j, r, c, pr, pc, rank, remr, sz;
    slong nr = A->r;
    slong nc = A->c;
    gr_lil_mat_t Rt;
    gr_sparse_vec_struct *pcol, *prow, *row, *col;
    gr_ptr cinv, cc;
    int status = GR_SUCCESS;

    if (gr_mat_is_compatible(A, R, ctx) != T_TRUE)
    {
        return GR_DOMAIN;
    }
    if (nr == 0 || nc == 0)
    {
        *res_rank = 0;
        return GR_SUCCESS;
    }
    if (A->nnz == 0)
    {
        gr_lil_mat_zero(R, ctx);
        return GR_SUCCESS;
    }
    // Get transpose and copy of A
    gr_lil_mat_init(Rt, nc, nr, ctx);
    status |= gr_lil_mat_transpose(Rt, A, ctx);
    status |= gr_lil_mat_set(R, A, ctx);

    GR_TMP_INIT2(cinv, cc, ctx);

    sz = ctx->sizeof_elem;

    /* Set up permutations */
    P = flint_malloc(nr*sizeof(*P));
    remr = nr;
    for (r = 0; r < nr; ++r) 
    {
        if (!A->rows[r].nnz || A->rows[r].inds[0] >= nc) P[r] = --remr; 
        else P[r] = -1;
    }
    
    /* Run elimination */
    rank = 0;
    for (pc = 0; pc < nc; ++pc)
    {
        pcol = &Rt->rows[pc];

        /* Get lowest weight incident row not used as previous pivot */
        pr = -1, prow = NULL;
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->inds[j], row = &R->rows[r];
            if (P[r] >= 0) continue;
            if (pr==-1 || (row->nnz < prow->nnz)) pr = r, prow = row;
        }
        if (pr == -1) continue;
        P[pr] = rank; 

        status |= gr_inv(cinv, gr_sparse_vec_find_entry(prow, pc, ctx), ctx);
        status |= gr_sparse_vec_mul_scalar(prow, prow, cinv, ctx);

        /* Gaussian eliminate rows */
        for (j = 0; j < pcol->nnz; ++j)
        {
            r = pcol->inds[j], row = &R->rows[r];
            if (r==pr) {status |= gr_zero(GR_ENTRY(pcol->nzs, j, sz), ctx); continue;}

            status |= gr_neg(cc, gr_sparse_vec_find_entry(row, pc, ctx), ctx);
            status |= gr_sparse_vec_addmul_scalar(row, prow, cc, ctx);
            if (row->nnz == 0 || row->inds[0] >= nc) P[r] = --remr;
        }
        /* Gaussian eliminate cols */
        status |= gr_sparse_vec_mul_scalar(pcol, pcol, cinv, ctx);
        for (j = 0; j < prow->nnz; ++j)
        {
            c = prow->inds[j], col = &Rt->rows[c];
            if (c >= nc || c==pc) continue;
            status |= gr_neg(cc, gr_sparse_vec_find_entry(col, pr, ctx), ctx);
            status |= gr_sparse_vec_addmul_scalar(col, pcol, cc, ctx);
        }
        rank += 1;
    }
    gr_lil_mat_clear(Rt, ctx);

    /* Fix nnz */
    R->nnz = 0;
    for (j = 0; j < nr; ++j)
    {
        R->nnz += R->rows[j].nnz;
    }

    /* Reorder rows */
    status |= gr_lil_mat_permute_rows(R, P, ctx);
    flint_free(P);
    GR_TMP_CLEAR2(cinv, cc, ctx);
    *res_rank = rank;
    return status;
}
