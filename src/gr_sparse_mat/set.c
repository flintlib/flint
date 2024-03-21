/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

typedef struct
{
    slong i;
    slong row;
    slong col;
    void * entry;
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

sparse_mat_index_t * _sort_coords(ulong * rows, ulong * cols, const void * entries, slong sz, slong num)
{
    slong i;
    sparse_mat_index_t * si;

    si = flint_malloc(num * sizeof(sparse_mat_index_t));
    for (i = 0; i < num; i++)
    {
        si[i].i = i;
        si[i].row = rows[i];
        si[i].col = cols[i];
        si[i].entry = (void *) (((char *) (entries)) + ((i) * (sz)));
    }

    qsort(si, num, sizeof(sparse_mat_index_t), sparse_mat_index_cmp);
    return si;
}

int gr_csr_mat_set(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx) {
    int status = GR_SUCCESS;

    if (dst != src)
    {
        if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        {
            return GR_DOMAIN;
        }
        dst->nnz = src->nnz;
        gr_csr_mat_fit_nnz(dst, src->nnz, ctx);
        memcpy(dst->rows, src->rows, src->r * sizeof(ulong));
        memcpy(dst->cols, src->cols, src->nnz * sizeof(ulong));
        status = _gr_vec_set(dst->nzs, src->nzs, src->nnz, ctx);
    }
    return status;
}

int gr_lil_mat_set(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx) {
    ulong row;
    int status = GR_SUCCESS;

    if (dst != src)
    {
        if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        {
            return GR_DOMAIN;
        }
        dst->nnz = src->nnz;

        for (row = 0; row < src->r; ++row) {
            status |= gr_sparse_vec_set(&dst->rows[row], &src->rows[row], ctx);
        }
    }
    return status;
}

int gr_coo_mat_set(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx) {
    int status = GR_SUCCESS;

    if (dst != src)
    {
        if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        {
            return GR_DOMAIN;
        }
        dst->nnz = src->nnz;
        gr_coo_mat_fit_nnz(dst, src->nnz, ctx);
        memcpy(dst->rows, src->rows, src->nnz * sizeof(ulong));
        memcpy(dst->cols, src->cols, src->nnz * sizeof(ulong));
        status = _gr_vec_set(dst->nzs, src->nzs, src->nnz, ctx);
        dst->is_canonical = src->is_canonical;
    }
    return status;
}

int gr_csr_mat_set_lil_mat(gr_csr_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    ulong row;
    int status = GR_SUCCESS;
    gr_sparse_vec_t dst_row;

    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    gr_csr_mat_fit_nnz(dst, src->nnz, ctx);

    dst->rows[0] = 0;
    for(row = 0; row < src->r; row++) {
        dst->rows[row+1] = dst->rows[row] + src->rows[row].nnz;
        _gr_csr_mat_borrow_row(dst_row, dst, row, ctx);
        status |= gr_sparse_vec_set(dst_row, &src->rows[row], ctx);
    }
    dst->nnz = src->nnz;

    return status;
}

int gr_lil_mat_set_csr_mat(gr_lil_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx)
{
    ulong row;
    int status = GR_SUCCESS;
    gr_sparse_vec_t mat_row;
    
    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    for(row = 0; row < src->r; row++) {
        _gr_csr_mat_borrow_row(mat_row, src, row, ctx);
        status |= gr_sparse_vec_set(&dst->rows[row], mat_row, ctx);
    }
    dst->nnz = src->nnz;

    return status;
}

int gr_csr_mat_set_coo_mat(gr_csr_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx)
{
   // Sort entries by row and column
    slong i, sz, nnz, row;
    int status = GR_SUCCESS;
    sparse_mat_index_t * si;
    gr_ptr entry;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);

    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    if (src->nnz == 0)
    {
        gr_csr_mat_zero(dst, ctx);
        return GR_SUCCESS;
    }
    
    sz = ctx->sizeof_elem;

    gr_csr_mat_fit_nnz(dst, src->nnz, ctx);
    if (src->is_canonical == T_TRUE)
    {
        // Just copy over the data and set the row offsets
        dst->rows[0] = 0;
        row = 0;
        for (i = 0; i < src->nnz; ++i)
            while (row != src->rows[i]) 
                dst->rows[++row] = i;
        while (row < dst->r)
            dst->rows[++row] = i;

        memcpy(dst->cols, src->cols, src->nnz * sizeof(ulong));
        status |= _gr_vec_set(dst->nzs, src->nzs, src->nnz, ctx);
        //_gr_vec_print(dst->nzs, src->nnz, ctx);
        dst->nnz = src->nnz;
    }
    else
    {
        // Sort coordinates
        si = _sort_coords(src->rows, src->cols, src->nzs, sz, src->nnz);
        
        // Accumulate nonzeroes into matrix
        row = 0;
        nnz = 0;
        dst->rows[0] = 0;
        entry = NULL;
        for (i = 0; i < src->nnz; ++i)
        {   
            //flint_printf("(%d, %d)\n", si[i].row, si[i].col);

            // Check if we can just accumulate
            if(i > 0 && si[i-1].row == si[i].row && si[i-1].col == si[i].col)
                status |= gr_add(entry, entry, si[i].entry, ctx);
            else
            {
                // Advance row offsets as needed

                // If previous entry does not exist or is not zero, advance to the next one
                if (entry == NULL || is_zero(entry, ctx) != T_TRUE)
                {
                    entry = GR_ENTRY(dst->nzs, nnz, sz);
                    ++nnz;
                }
                while (row != si[i].row) 
                    dst->rows[++row] = nnz - 1;
                dst->cols[nnz - 1] = si[i].col;
                status |= gr_set(entry, si[i].entry, ctx);
            }
        }
        if (entry != NULL && is_zero(entry, ctx) == T_TRUE)
            --nnz;

        // Set remaining row offsets and overall number of nonzeroes
        while (row < dst->r)
            dst->rows[++row] = nnz;
        dst->nnz = nnz;
    }
    return status;
}

int gr_lil_mat_set_coo_mat(gr_lil_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx)
{
    slong i, sz, nnz, row, row_start_idx, row_end_idx;
    int status = GR_SUCCESS;
    sparse_mat_index_t * si;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);
    gr_ptr entry;
    
    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    if (src->nnz == 0)
    {
        gr_lil_mat_zero(dst, ctx);
        return GR_SUCCESS;
    }

    sz = ctx->sizeof_elem;
    if (src->is_canonical == T_TRUE)
    {
        row_start_idx = row_end_idx = 0;
        for (row = 0; row < dst->r; ++row)
        {
            while (row_end_idx < src->nnz && src->rows[row_end_idx] == row) ++row_end_idx;
            status |= gr_sparse_vec_from_entries(
                &dst->rows[row], src->cols + row_start_idx, GR_ENTRY(src->nzs, row_start_idx, sz), row_end_idx - row_start_idx, T_TRUE, ctx
            );
            row_start_idx = row_end_idx;
        }
        dst->nnz = src->nnz;
    }
    else
    {
        // Sort entries by row and column
        si = _sort_coords(src->rows, src->cols, src->nzs, sz, src->nnz);

        // Construct rows one by one
        row_start_idx = 0;
        dst->nnz = 0;
        for (row = 0; row < dst->r; ++row)
        {
            // Get range of indicies in current row and estimate of nonzeroes
            for (row_end_idx = row_start_idx; row_end_idx < src->nnz; ++row_end_idx)
                if (si[row_end_idx].row != row)
                    break;
            //flint_printf("(%d, %d)\n", row_start_idx, row_end_idx);
            gr_sparse_vec_fit_nnz(&dst->rows[row], row_end_idx - row_start_idx, ctx);

            // Add nonzeroes to row
            entry = NULL;
            nnz = 0;
            for (i = row_start_idx; i < row_end_idx; ++i)
            {
                //flint_printf("\t(%d, %d)\n", si[i].row, si[i].col);
                // Skip zero entries
                if (is_zero(si[i].entry, ctx) == T_TRUE)
                    continue;

                // Check if we can just accumulate
                if(i > row_start_idx && si[i-1].col == si[i].col)
                    status |= gr_add(entry, entry, si[i].entry, ctx);
                else 
                {
                    // Check if need to get new entry to store to
                    if (entry == NULL || is_zero(entry, ctx) != T_TRUE)
                    {
                        entry = GR_ENTRY(dst->rows[row].nzs, nnz, sz);
                        ++nnz;
                    }
                    // Store to entry and update current column
                    dst->rows[row].inds[nnz - 1] = si[i].col;
                    status |= gr_set(entry, si[i].entry, ctx);
                }
            }
            if (entry != NULL && is_zero(entry, ctx) == T_TRUE)
                --nnz;
            dst->rows[row].nnz = nnz;
            dst->nnz += nnz;
            row_start_idx = row_end_idx;
        }
    }

    return status;
}

int gr_coo_mat_set_lil_mat(gr_coo_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    ulong row, i;
    slong sz;
    int status = GR_SUCCESS;
    gr_sparse_vec_struct *vec;

    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    gr_coo_mat_fit_nnz(dst, src->nnz, ctx);

    sz = ctx->sizeof_elem;

    dst->nnz = 0;
    for(row = 0; row < src->r; row++) {
        vec = &src->rows[row];
        for (i = 0; i < vec->nnz; ++i)
            dst->rows[dst->nnz + i] = row;
        
        memcpy(dst->cols + dst->nnz, vec->inds, vec->nnz * sizeof(ulong));
        status |= _gr_vec_set(GR_ENTRY(dst->nzs, dst->nnz, sz), vec->nzs, vec->nnz, ctx);
        dst->nnz += vec->nnz;
    }
    dst->nnz = src->nnz;
    dst->is_canonical = T_TRUE;
    return status;
}

int gr_coo_mat_set_csr_mat(gr_coo_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx)
{
    ulong row;
    int status = GR_SUCCESS;
    
    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    gr_coo_mat_fit_nnz(dst, src->nnz, ctx);

    dst->nnz = 0;
    for(row = 0; row < dst->r; ++row)
        while (dst->nnz < src->rows[row + 1])
            dst->rows[dst->nnz++] = row;
    memcpy(dst->cols, src->cols, src->nnz * sizeof(ulong));
    status |= _gr_vec_set(dst->nzs, src->nzs, src->nnz, ctx);
    dst->is_canonical = T_TRUE;

    return status;
}

int gr_csr_mat_set_mat(gr_csr_mat_t dst, const gr_mat_t src, gr_ctx_t ctx)
{
    ulong row, col;
    slong nnz, sz;
    int status = GR_SUCCESS;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);
    
    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    // Get number of nonzero entries
    sz = ctx->sizeof_elem;
    nnz = 0;
    for (row = 0; row < src->r; ++row)
        for (col = 0; col < src->c; ++col)
            if (is_zero(GR_MAT_ENTRY(src, row, col, sz), ctx) != T_TRUE)
                nnz++;
    gr_csr_mat_fit_nnz(dst, nnz, ctx);

    // Construct sparse matrix from nonzeroes
    dst->rows[0] = 0;
    nnz = 0;
    for(row = 0; row < src->r; row++)
    {
        for (col = 0; col < src->c; ++col)
        {
            if (is_zero(GR_MAT_ENTRY(src, row, col, sz), ctx) != T_TRUE)
            {
                dst->cols[nnz] = col;
                status |= gr_set(GR_ENTRY(dst->nzs, nnz, sz), GR_MAT_ENTRY(src, row, col, sz), ctx); 
                ++nnz;       
            }
        }
        dst->rows[row + 1] = nnz;
    }
    dst->nnz = nnz;

    return status;
}

int gr_lil_mat_set_mat(gr_lil_mat_t dst, const gr_mat_t src, gr_ctx_t ctx)
{
    ulong row;
    slong sz;
    int status = GR_SUCCESS;
    
    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    sz = ctx->sizeof_elem;
    dst->nnz = 0;
    for(row = 0; row < src->r; row++) {
        status |= gr_sparse_vec_set_vec(&dst->rows[row], GR_MAT_ENTRY(src, row, 0, sz), src->c, ctx);
        dst->nnz += dst->rows[row].nnz;
    }

    return status;
}

int gr_coo_mat_set_mat(gr_coo_mat_t dst, const gr_mat_t src, gr_ctx_t ctx)
{
    ulong row, col;
    slong sz, nnz;
    int status = GR_SUCCESS;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);
    
    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    sz = ctx->sizeof_elem;

    // Get number of nonzero entries
    nnz = 0;
    for (row = 0; row < src->r; ++row)
        for (col = 0; col < src->c; ++col)
            if (is_zero(GR_MAT_ENTRY(src, row, col, sz), ctx) != T_TRUE)
                nnz++;
    gr_coo_mat_fit_nnz(dst, nnz, ctx);

    // Construct sparse matrix from nonzeroes
    nnz = 0;
    for(row = 0; row < src->r; row++)
    {
        for (col = 0; col < src->c; ++col)
        {
            if (is_zero(GR_MAT_ENTRY(src, row, col, sz), ctx) != T_TRUE)
            {
                dst->rows[nnz] = row;
                dst->cols[nnz] = col;
                status |= gr_set(GR_ENTRY(dst->nzs, nnz, sz), GR_MAT_ENTRY(src, row, col, sz), ctx); 
                ++nnz;
            }
        }
    }
    dst->nnz = nnz;

    return status;
}
