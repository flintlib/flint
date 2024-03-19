/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "gr_sparse_mat.h"

typedef struct
{
    slong i;
    slong row;
    slong col;
}
sparse_mat_index_t;

static int sparse_mat_index_cmp(const void* a, const void* b)
{
    slong arow = ((sparse_mat_index_t*)(a))->row;
    slong brow = ((sparse_mat_index_t*)(b))->row;
    slong acol = ((sparse_mat_index_t*)(a))->col;
    slong bcol = ((sparse_mat_index_t*)(b))->col;
    return (arow < brow ? -1 : (arow > brow ? 1 : (acol < bcol ? -1 : (acol > bcol ? 1 : 0))));
}

static sparse_mat_index_t * _sort_coords(ulong * rows, ulong * cols, slong num)
{
    slong i;
    sparse_mat_index_t * si;

    si = flint_malloc(num * sizeof(sparse_mat_index_t));
    for (i = 0; i < num; i++)
    {
        si[i].i = i;
        si[i].row = rows[i];
        si[i].col = cols[i];
    }

    qsort(si, num, sizeof(sparse_mat_index_t), sparse_mat_index_cmp);
    return si;
}

int
gr_coo_mat_canonicalize(gr_coo_mat_t mat, gr_ctx_t ctx)
{
    slong i, j, k, sz, nnz;
    int status = GR_SUCCESS;
    sparse_mat_index_t *si;
    ulong *inv_si;
    gr_ptr entry;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);

    if (mat->is_canonical)
        return GR_SUCCESS;

    sz = ctx->sizeof_elem;
    
    // Get sorted order for matrices indices (and inverse mapping)
    si = _sort_coords(mat->rows, mat->cols, mat->nnz);
    inv_si = flint_malloc(mat->nnz * sizeof(ulong));
    for (i = 0; i < mat->nnz; ++i)
        inv_si[si[i].i] = i;

    // Use swaps to apply sort to entries
    for (i = 0; i < mat->nnz; ++i)
    {
        j = si[i].i;
        if (i != j)
        {
            FLINT_SWAP(ulong, mat->rows[i], mat->rows[j]);
            FLINT_SWAP(ulong, mat->cols[i], mat->cols[j]);
            gr_swap(GR_ENTRY(mat->nzs, i, sz), GR_ENTRY(mat->nzs, j, sz), ctx);

            // Fix mappings to remove i from permutation
            k = inv_si[i];
            si[k].i = j;
            inv_si[j] = k;
        }
    }
    flint_free(si);
    flint_free(inv_si);

    // Compress duplicated entries
    nnz = 0;
    entry = NULL;
    for (i = 0; i < mat->nnz; ++i)
    {
        if(i > 0 && mat->rows[i-1] == mat->rows[i] && mat->cols[i-1] == mat->cols[i])
            status |= gr_add(entry, entry, GR_ENTRY(mat->nzs, i, sz), ctx);
        else
        {
            // If previous entry does not exist or is not zero, advance to the next one
            if (entry == NULL || is_zero(entry, ctx) != T_TRUE)
            {
                entry = GR_ENTRY(mat->nzs, nnz, sz);
                ++nnz;
            }
            status |= gr_set(entry, GR_ENTRY(mat->nzs, i, sz), ctx);
        }
    }
    return status;
}


