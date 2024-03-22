/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef GR_SPARSE_MAT_H
#define GR_SPARSE_MAT_H

#ifdef GR_SPARSE_MAT_INLINES_C
#define GR_SPARSE_MAT_INLINE
#else
#define GR_SPARSE_MAT_INLINE static inline
#endif

#include "gr.h"
#include "gr_mat.h"
#include "gr_sparse_vec.h"
#include "gr_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

/**
 * Types and basic access
**/

typedef struct
{
    slong r;
    slong c;
    slong nnz;
    slong alloc;
    ulong * rows;
    ulong * cols;
    gr_ptr nzs;
}
gr_csr_mat_struct;

typedef gr_csr_mat_struct gr_csr_mat_t[1];

typedef struct
{
    slong r;
    slong c;
    slong nnz;
    gr_sparse_vec_struct * rows;    
}
gr_lil_mat_struct;

typedef gr_lil_mat_struct gr_lil_mat_t[1];

typedef struct
{
    slong r;
    slong c;
    slong nnz;
    slong alloc;
    ulong * rows;
    ulong * cols;
    gr_ptr nzs;
    truth_t is_canonical;
}
gr_coo_mat_struct;

typedef gr_coo_mat_struct gr_coo_mat_t[1];

#define gr_sparse_mat_nrows(mat, ctx) ((mat)->r)
#define gr_sparse_mat_ncols(mat, ctx) ((mat)->c)
#define gr_sparse_mat_nnz(mat, ctx) ((mat)->nnz)

#define GR_CSR_MAT_COL(mat,row,nz_idx) ((mat)->cols + (mat)->rows[row] + nz_idx)
#define GR_CSR_MAT_ENTRY(mat,row,nz_idx,sz) GR_ENTRY((mat)->nzs, (mat)->rows[row] + nz_idx, sz)

#define GR_LIL_MAT_COL(mat,row,nz_idx) ((mat)->rows[row].inds + nz_idx)
#define GR_LIL_MAT_ENTRY(mat,row,nz_idx,sz) GR_ENTRY((mat)->rows[row].nzs, nz_idx, sz)

#define GR_COO_MAT_ROW(mat,nz_idx) ((mat)->rows + nz_idx)
#define GR_COO_MAT_COL(mat,nz_idx) ((mat)->cols + nz_idx)
#define GR_COO_MAT_ENTRY(mat,nz_idx,sz) GR_ENTRY((mat)->nzs, nz_idx, sz)

GR_SPARSE_MAT_INLINE ulong * 
gr_csr_mat_col_ptr(gr_csr_mat_t mat, slong row, slong nz_idx)
{
    return GR_CSR_MAT_COL(mat, row, nz_idx);
}

GR_SPARSE_MAT_INLINE const ulong * 
gr_csr_mat_col_srcptr(const gr_csr_mat_t mat, slong row, slong nz_idx)
{
    if (row < 0 || row >= mat->r || nz_idx < 0 || mat->rows[row] + nz_idx >= mat->rows[row+1])
        return NULL;
    return GR_CSR_MAT_COL(mat, row, nz_idx);
}

GR_SPARSE_MAT_INLINE gr_ptr 
gr_csr_mat_entry_ptr(gr_csr_mat_t mat, slong row, slong nz_idx, gr_ctx_t ctx)
{
    if (row < 0 || row >= mat->r || nz_idx < 0 || mat->rows[row] + nz_idx >= mat->rows[row+1])
        return NULL;
    return GR_CSR_MAT_ENTRY(mat, row, nz_idx, ctx->sizeof_elem);
}

GR_SPARSE_MAT_INLINE gr_srcptr 
gr_csr_mat_entry_srcptr(const gr_csr_mat_t mat, slong row, slong nz_idx, gr_ctx_t ctx)
{
    if (row < 0 || row >= mat->r || nz_idx < 0 || mat->rows[row] + nz_idx >= mat->rows[row+1])
        return NULL;
    return GR_CSR_MAT_ENTRY(mat, row, nz_idx, ctx->sizeof_elem);
}

GR_SPARSE_MAT_INLINE ulong * 
gr_lil_mat_col_ptr(gr_lil_mat_t mat, slong row, slong nz_idx)
{
    if (row < 0 || row >= mat->r || nz_idx < 0 || nz_idx >= mat->rows[row].nnz)
        return NULL;
    return GR_LIL_MAT_COL(mat, row, nz_idx);
}

GR_SPARSE_MAT_INLINE const ulong * 
gr_lil_mat_col_srcptr(const gr_lil_mat_t mat, slong row, slong nz_idx)
{
    if (row < 0 || row >= mat->r || nz_idx < 0 || nz_idx >= mat->rows[row].nnz)
        return NULL;
    return GR_LIL_MAT_COL(mat, row, nz_idx);
}

GR_SPARSE_MAT_INLINE gr_ptr 
gr_lil_mat_entry_ptr(gr_lil_mat_t mat, slong row, slong nz_idx, gr_ctx_t ctx)
{
    if (row < 0 || row >= mat->r || nz_idx < 0 || nz_idx >= mat->rows[row].nnz)
        return NULL;
    return GR_LIL_MAT_ENTRY(mat, row, nz_idx, ctx->sizeof_elem);
}

GR_SPARSE_MAT_INLINE gr_srcptr 
gr_lil_mat_entry_srcptr(const gr_lil_mat_t mat, slong row, slong nz_idx, gr_ctx_t ctx)
{
    if (row < 0 || row >= mat->r || nz_idx < 0 || nz_idx >= mat->rows[row].nnz)
        return NULL;
    return GR_LIL_MAT_ENTRY(mat, row, nz_idx, ctx->sizeof_elem);
}

GR_SPARSE_MAT_INLINE ulong * 
gr_coo_mat_row_ptr(gr_coo_mat_t mat, slong nz_idx)
{
    if (nz_idx < 0 || nz_idx >= mat->nnz)
        return NULL;
    return GR_COO_MAT_ROW(mat, nz_idx);
}

GR_SPARSE_MAT_INLINE const ulong * 
gr_coo_mat_row_srcptr(const gr_coo_mat_t mat, slong nz_idx)
{
    if (nz_idx < 0 || nz_idx >= mat->nnz)
        return NULL;
    return GR_COO_MAT_ROW(mat, nz_idx);
}

GR_SPARSE_MAT_INLINE const ulong * 
gr_coo_mat_col_ptr(gr_coo_mat_t mat, slong nz_idx)
{
    if (nz_idx < 0 || nz_idx >= mat->nnz)
        return NULL;
    return GR_COO_MAT_COL(mat, nz_idx);
}

GR_SPARSE_MAT_INLINE const ulong * 
gr_coo_mat_col_srcptr(const gr_coo_mat_t mat, slong nz_idx)
{
    if (nz_idx < 0 || nz_idx >= mat->nnz)
        return NULL;
    return GR_COO_MAT_COL(mat, nz_idx);
}

GR_SPARSE_MAT_INLINE gr_ptr 
gr_coo_mat_entry_ptr(gr_coo_mat_t mat, slong nz_idx, gr_ctx_t ctx)
{
    if (nz_idx < 0 || nz_idx >= mat->nnz)
        return NULL;
    return GR_COO_MAT_ENTRY(mat, nz_idx, ctx->sizeof_elem);
}

GR_SPARSE_MAT_INLINE gr_srcptr
gr_coo_mat_entry_srcptr(const gr_coo_mat_t mat, slong nz_idx, gr_ctx_t ctx)
{
    if (nz_idx < 0 || nz_idx >= mat->nnz)
        return NULL;
    return GR_COO_MAT_ENTRY(mat, nz_idx, ctx->sizeof_elem);
}

/* Generics */
/*
typedef int ((*gr_method_mat_unary_op_get_scalar)(gr_ptr, const gr_sparse_mat_t, gr_ctx_ptr));
typedef int ((*gr_method_mat_unary_op)(gr_sparse_mat_t, const gr_sparse_mat_t, gr_ctx_ptr));
typedef int ((*gr_method_mat_binary_op)(gr_sparse_mat_t, const gr_sparse_mat_t, const gr_sparse_mat_t, gr_ctx_ptr));
typedef int ((*gr_method_mat_pivot_op)(slong *, gr_sparse_mat_t, slong, slong, slong, gr_ctx_ptr));
typedef int ((*gr_method_mat_diagonalization_op)(gr_vec_t, gr_sparse_mat_t, gr_sparse_mat_t, const gr_sparse_mat_t, int, gr_ctx_ptr));

#define GR_SPARSE_MAT_UNARY_OP_GET_SCALAR(ctx, NAME) (((gr_method_mat_unary_op_get_scalar *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_MAT_UNARY_OP(ctx, NAME) (((gr_method_mat_unary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_MAT_BINARY_OP(ctx, NAME) (((gr_method_mat_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_MAT_PIVOT_OP(ctx, NAME) (((gr_method_mat_pivot_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_MAT_DIAGONALIZATION_OP(ctx, NAME) (((gr_method_mat_diagonalization_op *) ctx->methods)[GR_METHOD_ ## NAME])
*/

GR_SPARSE_MAT_INLINE void
gr_csr_mat_init(gr_csr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx) {
    memset(mat, 0, sizeof(gr_csr_mat_t));
    mat->r = rows;
    mat->c = cols;
    mat->rows = flint_calloc(rows + 1, sizeof(ulong));
}

GR_SPARSE_MAT_INLINE void
gr_lil_mat_init(gr_lil_mat_t mat, slong rows, slong cols, gr_ctx_t ctx) {
    int row;

    memset(mat, 0, sizeof(gr_lil_mat_t));
    mat->r = rows;
    mat->c = cols;
    mat->rows = flint_calloc(rows, sizeof(gr_sparse_vec_struct));

    for(row = 0; row < mat->r; ++row)
        gr_sparse_vec_init(&mat->rows[row], cols, ctx);
}

GR_SPARSE_MAT_INLINE void
gr_coo_mat_init(gr_coo_mat_t mat, slong rows, slong cols, gr_ctx_t ctx) {
    memset(mat, 0, sizeof(gr_coo_mat_t));
    mat->r = rows;
    mat->c = cols;
    mat->is_canonical = T_TRUE;
}

GR_SPARSE_MAT_INLINE void
gr_csr_mat_clear(gr_csr_mat_t mat, gr_ctx_t ctx) {
    if (mat->alloc != 0)
    {
        _gr_vec_clear(mat->nzs, mat->alloc, ctx);
        flint_free(mat->nzs);
        flint_free(mat->cols);
    }
    flint_free(mat->rows);
    memset(mat, 0, sizeof(gr_csr_mat_t));
}

GR_SPARSE_MAT_INLINE void
gr_lil_mat_clear(gr_lil_mat_t mat, gr_ctx_t ctx) {
    int row;

    for (row = 0; row < mat->r; ++row)
        gr_sparse_vec_clear(&mat->rows[row], ctx);
    flint_free(mat->rows);
    memset(mat, 0, sizeof(gr_lil_mat_t));
}

GR_SPARSE_MAT_INLINE void
gr_coo_mat_clear(gr_coo_mat_t mat, gr_ctx_t ctx) {
    if (mat->alloc != 0)
    {
        _gr_vec_clear(mat->nzs, mat->alloc, ctx);
        flint_free(mat->nzs);
        flint_free(mat->rows);
        flint_free(mat->cols);
    }
    memset(mat, 0, sizeof(gr_csr_mat_t));
}

GR_SPARSE_MAT_INLINE void
gr_csr_mat_swap(gr_csr_mat_t mat1, gr_csr_mat_t mat2, gr_ctx_t ctx)
{
    FLINT_SWAP(gr_csr_mat_struct, *mat1, *mat2);
}

GR_SPARSE_MAT_INLINE void
gr_lil_mat_swap(gr_lil_mat_t mat1, gr_lil_mat_t mat2, gr_ctx_t ctx)
{
    FLINT_SWAP(gr_lil_mat_struct, *mat1, *mat2);
}

GR_SPARSE_MAT_INLINE void
gr_coo_mat_swap(gr_coo_mat_t mat1, gr_coo_mat_t mat2, gr_ctx_t ctx)
{
    FLINT_SWAP(gr_coo_mat_struct, *mat1, *mat2);
}

void gr_csr_mat_fit_nnz(gr_csr_mat_t mat, slong nnz, gr_ctx_t ctx);
void gr_lil_mat_fit_nnz(gr_coo_mat_t mat, slong *nnz, gr_ctx_t ctx);
void gr_coo_mat_fit_nnz(gr_coo_mat_t mat, slong nnz, gr_ctx_t ctx);


void gr_csr_mat_shrink_to_nnz(gr_csr_mat_t mat, gr_ctx_t ctx);
void gr_lil_mat_shrink_to_nnz(gr_lil_mat_t mat, gr_ctx_t ctx);
void gr_coo_mat_shrink_to_nnz(gr_coo_mat_t mat, gr_ctx_t ctx);


void gr_csr_mat_set_cols(gr_csr_mat_t mat, slong cols, gr_ctx_t ctx);
void gr_lil_mat_set_cols(gr_lil_mat_t mat, slong cols, gr_ctx_t ctx);
void gr_coo_mat_set_cols(gr_coo_mat_t mat, slong cols, gr_ctx_t ctx);

int gr_coo_mat_from_entries(gr_coo_mat_t mat, ulong *rows, ulong *cols, gr_srcptr entries, slong nnz, truth_t is_canonical, gr_ctx_t ctx);

truth_t gr_coo_mat_is_canonical(gr_coo_mat_t mat, gr_ctx_t ctx);

int gr_coo_mat_canonicalize(gr_coo_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int 
gr_coo_mat_randtest(gr_coo_mat_t mat, slong nnz, int replacement, truth_t is_canonical, flint_rand_t state, gr_ctx_t ctx);

WARN_UNUSED_RESULT int 
gr_coo_mat_randtest_prob(gr_coo_mat_t mat, double prob, flint_rand_t state, gr_ctx_t ctx);

/**
 * Getting, setting, and conversion
**/

GR_SPARSE_MAT_INLINE void
_gr_csr_mat_borrow_row(gr_sparse_vec_t res, const gr_csr_mat_t mat, slong r, gr_ctx_t ctx) {
    ulong offset;

    offset = mat->rows[r];
    res->length = mat->c;
    res->nnz = mat->rows[r+1] - offset;
    res->alloc = res->nnz;
    res->inds = mat->cols + offset;
    res->nzs = GR_ENTRY(mat->nzs, offset, ctx->sizeof_elem);
}

GR_SPARSE_MAT_INLINE void
_gr_coo_mat_borrow_row(gr_sparse_vec_t res, gr_coo_mat_t mat, slong r, gr_ctx_t ctx) {
    ulong offset, nnz;

    if (mat->is_canonical == T_FALSE)
        gr_coo_mat_canonicalize(mat, ctx);
    
    // Binary search to find start and end of row
    offset = 0;
    nnz = 0;

    res->length = mat->c;
    res->nnz = nnz;
    res->alloc = nnz;
    res->inds = mat->cols + offset;
    res->nzs = GR_ENTRY(mat->nzs, offset, ctx->sizeof_elem);
}

gr_ptr gr_csr_mat_find_entry(gr_csr_mat_t mat, slong row, slong col, gr_ctx_t ctx);
gr_ptr gr_lil_mat_find_entry(gr_lil_mat_t mat, slong row, slong col, gr_ctx_t ctx);
gr_ptr gr_coo_mat_find_entry(gr_coo_mat_t mat, slong row, slong col, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_csr_mat_get_entry(gr_ptr res, gr_csr_mat_t mat, slong row, slong col, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_get_entry(gr_ptr res, gr_lil_mat_t mat, slong row, slong col, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_get_entry(gr_ptr res, gr_coo_mat_t mat, slong row, slong col, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_csr_mat_set_entry(gr_csr_mat_t mat, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_set_entry(gr_lil_mat_t mat, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_set_entry(gr_coo_mat_t mat, slong row, slong col, gr_srcptr entry, gr_ctx_t ctx);


GR_SPARSE_MAT_INLINE void 
gr_csr_mat_zero(gr_csr_mat_t mat, gr_ctx_t ctx) {
    mat->nnz = 0;
    memset(mat->rows, 0, (mat->r + 1) * sizeof(ulong));
}

GR_SPARSE_MAT_INLINE void 
gr_lil_mat_zero(gr_lil_mat_t mat, gr_ctx_t ctx) {
    int row;

    mat->nnz = 0;
    for(row = 0; row < mat->r; ++row)
    {
        gr_sparse_vec_zero(&mat->rows[row], ctx);
    }
}

GR_SPARSE_MAT_INLINE void 
gr_coo_mat_zero(gr_coo_mat_t mat, gr_ctx_t ctx) {
    mat->nnz = 0;
}

WARN_UNUSED_RESULT int
gr_csr_mat_set(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int 
gr_lil_mat_set(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int 
gr_coo_mat_set(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_csr_mat_set_lil_mat(gr_csr_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_lil_mat_set_csr_mat(gr_lil_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_coo_mat_set_csr_mat(gr_coo_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_coo_mat_set_lil_mat(gr_coo_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_lil_mat_set_coo_mat(gr_lil_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_csr_mat_set_coo_mat(gr_csr_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_csr_mat_set_mat(gr_csr_mat_t dst, const gr_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_lil_mat_set_mat(gr_lil_mat_t dst, const gr_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_coo_mat_set_mat(gr_coo_mat_t dst, const gr_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_mat_set_csr_mat(gr_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_mat_set_lil_mat(gr_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_mat_set_coo_mat(gr_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx);

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_csr_mat_init_set(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx) {
    gr_csr_mat_init(dst, src->r, src->c, ctx);
    return gr_csr_mat_set(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_lil_mat_init_set(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx) {
    gr_lil_mat_init(dst, src->r, src->c, ctx);
    return gr_lil_mat_set(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_coo_mat_init_set(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx) {
    gr_coo_mat_init(dst, src->r, src->c, ctx);
    return gr_coo_mat_set(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_csr_mat_init_set_lil_mat(gr_csr_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    gr_csr_mat_init(dst, src->r, src->c, ctx);
    return gr_csr_mat_set_lil_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_lil_mat_init_set_csr_mat(gr_lil_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx)
{
    gr_lil_mat_init(dst, src->r, src->c, ctx);
    return gr_lil_mat_set_csr_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_csr_mat_init_set_coo_mat(gr_csr_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx)
{
    gr_csr_mat_init(dst, src->r, src->c, ctx);
    return gr_csr_mat_set_coo_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_lil_mat_init_set_coo_mat(gr_lil_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx)
{
    gr_lil_mat_init(dst, src->r, src->c, ctx);
    return gr_lil_mat_set_coo_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_coo_mat_init_set_csr_mat(gr_coo_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx)
{
    gr_coo_mat_init(dst, src->r, src->c, ctx);
    return gr_coo_mat_set_csr_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_coo_mat_init_set_lil_mat(gr_coo_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    gr_coo_mat_init(dst, src->r, src->c, ctx);
    return gr_coo_mat_set_lil_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_csr_mat_init_set_mat(gr_csr_mat_t dst, const gr_mat_t src, gr_ctx_t ctx)
{
    gr_csr_mat_init(dst, src->r, src->c, ctx);
    return gr_csr_mat_set_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_lil_mat_init_set_mat(gr_lil_mat_t dst, const gr_mat_t src, gr_ctx_t ctx)
{
    gr_lil_mat_init(dst, src->r, src->c, ctx);
    return gr_lil_mat_set_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_coo_mat_init_set_mat(gr_coo_mat_t dst, const gr_mat_t src, gr_ctx_t ctx)
{
    gr_coo_mat_init(dst, src->r, src->c, ctx);
    return gr_coo_mat_set_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_mat_init_set_csr_mat(gr_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx)
{
    gr_mat_init(dst, src->r, src->c, ctx);
    return gr_mat_set_csr_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_mat_init_set_lil_mat(gr_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    gr_mat_init(dst, src->r, src->c, ctx);
    return gr_mat_set_lil_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int
gr_mat_init_set_coo_mat(gr_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx)
{
    gr_mat_init(dst, src->r, src->c, ctx);
    return gr_mat_set_coo_mat(dst, src, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int 
gr_lil_mat_update(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    slong row;
    int status = GR_SUCCESS;

    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
        return GR_DOMAIN;

    for (row = 0; row < dst->r; ++row)
        status |= gr_sparse_vec_update(&dst->rows[row], &src->rows[row], ctx);
    return status;
}

WARN_UNUSED_RESULT int gr_csr_mat_permute_cols(gr_csr_mat_t mat, slong * perm, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_permute_cols(gr_lil_mat_t mat, slong * perm, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_permute_cols(gr_coo_mat_t mat, slong * perm, gr_ctx_t ctx);

void gr_lil_mat_window_init(gr_lil_mat_t window, const gr_lil_mat_t mat, slong r1, slong c1, slong r2, slong c2, gr_ctx_t ctx);

GR_SPARSE_MAT_INLINE void
gr_lil_mat_window_clear(gr_lil_mat_t window, gr_ctx_t ctx)
{
    flint_free(window->rows);
}

WARN_UNUSED_RESULT int gr_lil_mat_swap_rows(gr_lil_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_swap_rows(gr_coo_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_lil_mat_permute_rows(gr_lil_mat_t mat, const slong * perm, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_permute_rows(gr_coo_mat_t mat, const slong * perm, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_csr_mat_invert_rows(gr_csr_mat_t mat, slong * perm, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_invert_rows(gr_lil_mat_t mat, slong * perm, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_invert_rows(gr_coo_mat_t mat, slong * perm, gr_ctx_t ctx);

/*
WARN_UNUSED_RESULT int gr_lil_mat_concat_horizontal(gr_lil_mat_t res, const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_concat_vertical(gr_csr_mat_t res, const gr_csr_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_concat_vertical(gr_lil_mat_t res, const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx);
*/

/*
WARN_UNUSED_RESULT int gr_sparse_mat_randops(gr_csr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_randpermdiag(int * parity, gr_csr_mat_t mat, flint_rand_t state, gr_ptr diag, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_randrank(gr_csr_mat_t mat, flint_rand_t state, slong rank, gr_ctx_t ctx);
*/

/**
 * Comparison
**/

GR_SPARSE_MAT_INLINE truth_t 
gr_csr_mat_is_zero(const gr_csr_mat_t mat, gr_ctx_t ctx)
{ return _gr_vec_is_zero(mat->nzs, mat->nnz, ctx); }

GR_SPARSE_MAT_INLINE truth_t 
gr_lil_mat_is_zero(const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    slong row;
    truth_t row_is_zero;
    truth_t ret = T_TRUE;

    for (row = 0; row < mat->r; ++row)
    {
        row_is_zero = gr_sparse_vec_is_zero(&mat->rows[row], ctx);
        if (row_is_zero == T_FALSE)
            return T_FALSE;
        else if (row_is_zero == T_UNKNOWN)
            ret = T_UNKNOWN;
    }
    return ret;
}

GR_SPARSE_MAT_INLINE truth_t 
gr_coo_mat_is_zero(gr_coo_mat_t mat, gr_ctx_t ctx)
{
    if (mat->is_canonical == T_FALSE)
        gr_coo_mat_canonicalize(mat, ctx);
    return _gr_vec_is_zero(mat->nzs, mat->nnz, ctx);
}

truth_t gr_csr_mat_is_one(const gr_csr_mat_t mat, gr_ctx_t ctx);
truth_t gr_lil_mat_is_one(const gr_lil_mat_t mat, gr_ctx_t ctx);
truth_t gr_coo_mat_is_one(const gr_lil_mat_t mat, gr_ctx_t ctx);

truth_t gr_csr_mat_is_neg_one(const gr_csr_mat_t mat, gr_ctx_t ctx);
truth_t gr_lil_mat_is_neg_one(const gr_lil_mat_t mat, gr_ctx_t ctx);
truth_t gr_coo_mat_is_neg_one(const gr_lil_mat_t mat, gr_ctx_t ctx);

truth_t gr_csr_mat_equal(const gr_csr_mat_t mat1, const gr_csr_mat_t mat2, gr_ctx_t ctx);
truth_t gr_lil_mat_equal(const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx);
truth_t gr_coo_mat_equal(const gr_coo_mat_t mat1, const gr_coo_mat_t mat2, gr_ctx_t ctx);
truth_t gr_csr_mat_equal_lil_mat(const gr_csr_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx);

/*
GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_csr_mat_one(gr_csr_mat_t res, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_set_scalar(gr_csr_mat_t res, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_set_ui(gr_csr_mat_t res, ulong v, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_set_si(gr_csr_mat_t res, slong v, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_set_fmpz(gr_csr_mat_t res, const fmpz_t v, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_set_fmpq(gr_csr_mat_t res, const fmpq_t v, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_csr_mat_set_fmpz_csr_mat(gr_csr_mat_t res, const fmpz_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_set_fmpq_csr_mat(gr_csr_mat_t res, const fmpq_csr_mat_t mat, gr_ctx_t ctx);
*/

/**
 * Output
**/

WARN_UNUSED_RESULT int gr_csr_mat_write_nz(gr_stream_t out, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_write_nz(gr_stream_t out, const gr_lil_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_write_nz(gr_stream_t out, const gr_coo_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_print_nz(const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_print_nz(const gr_lil_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_print_nz(const gr_coo_mat_t mat, gr_ctx_t ctx);

/**
 * Arithmetic
**/

WARN_UNUSED_RESULT int gr_csr_mat_neg(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_neg(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_neg(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_lil_mat_add(gr_lil_mat_t dst, const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_sub(gr_lil_mat_t dst, const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx);
//WARN_UNUSED_RESULT int gr_lil_mat_mul(gr_lil_mat_t dst, const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_lil_mat_addmul_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_submul_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx);

/**
 * Component-wise multiplication and division
**/

WARN_UNUSED_RESULT int gr_csr_mat_mul_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_mul_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_mul_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_mul_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_mul_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_mul_scalar_2exp_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_div_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_div_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_div_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_div_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_div_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_divexact_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_divexact_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_divexact_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_divexact_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_csr_mat_divexact_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_lil_mat_mul_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_mul_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_mul_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_mul_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_mul_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_mul_scalar_2exp_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_div_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_div_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_div_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_div_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_div_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_divexact_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_divexact_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_divexact_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_divexact_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_divexact_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_coo_mat_mul_scalar(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_mul_scalar_si(gr_coo_mat_t dst, const gr_coo_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_mul_scalar_ui(gr_coo_mat_t dst, const gr_coo_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_mul_scalar_fmpz(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_mul_scalar_fmpq(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_mul_scalar_2exp_si(gr_coo_mat_t dst, const gr_coo_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_div_scalar(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_div_scalar_si(gr_coo_mat_t dst, const gr_coo_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_div_scalar_ui(gr_coo_mat_t dst, const gr_coo_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_div_scalar_fmpz(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_div_scalar_fmpq(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_divexact_scalar(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_divexact_scalar_si(gr_coo_mat_t dst, const gr_coo_mat_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_divexact_scalar_ui(gr_coo_mat_t dst, const gr_coo_mat_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_divexact_scalar_fmpz(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_coo_mat_divexact_scalar_fmpq(gr_coo_mat_t dst, const gr_coo_mat_t src, const fmpq_t c, gr_ctx_t ctx);

/**
 * Sum and product
**/

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_csr_mat_sum(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx)
{ return _gr_vec_sum(res, mat->nzs, mat->nnz, ctx); }

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_lil_mat_sum(gr_ptr res, const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    slong row;
    int status = GR_SUCCESS;
    gr_ptr elem;

    elem = flint_malloc(ctx->sizeof_elem);
    gr_init(elem, ctx);
    for (row = 0; row < mat->r; ++row)
    {
        status |= gr_sparse_vec_sum(elem, &mat->rows[row], ctx);
        status |= gr_add(res, res, elem, ctx);
    }
    gr_clear(elem, ctx);
    return status;
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_coo_mat_sum(gr_ptr res, const gr_coo_mat_t mat, gr_ctx_t ctx)
{ return _gr_vec_sum(res, mat->nzs, mat->nnz, ctx); }

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_csr_mat_nz_product(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx)
{ return _gr_vec_product(res, mat->nzs, mat->nnz, ctx); }

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_lil_mat_nz_product(gr_ptr res, const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    slong row;
    int status = GR_SUCCESS;
    gr_ptr elem;

    elem = flint_malloc(ctx->sizeof_elem);
    gr_init(elem, ctx);
    status |= gr_one(elem, ctx);
    for (row = 0; row < mat->r; ++row)
    {
        status |= gr_sparse_vec_nz_product(elem, &mat->rows[row], ctx);
        status |= gr_mul(res, res, elem, ctx);
    }
    gr_clear(elem, ctx);
    return status;
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_coo_mat_nz_product(gr_ptr res, const gr_coo_mat_t mat, gr_ctx_t ctx)
{ return _gr_vec_product(res, mat->nzs, mat->nnz, ctx); }

/**
 * Transpose
*/
int gr_lil_mat_transpose(gr_lil_mat_t B, const gr_lil_mat_t A, gr_ctx_t ctx);

/**
 * Matrix multiplication
**/

WARN_UNUSED_RESULT int gr_csr_mat_mul_vec(gr_ptr v, const gr_csr_mat_t A, gr_srcptr u, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_mul_vec(gr_ptr v, const gr_lil_mat_t A, gr_srcptr u, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_csr_mat_mul_mat_transpose(gr_mat_t Ct, const gr_csr_mat_t A, const gr_mat_t Bt, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_mul_mat_transpose(gr_mat_t Ct, const gr_lil_mat_t A, const gr_mat_t Bt, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_csr_mat_mul_mat(gr_mat_t C, const gr_csr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_mul_mat(gr_mat_t C, const gr_lil_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

/**
 * Solving, nullvector, and nullspace computation
**/

WARN_UNUSED_RESULT int gr_lil_mat_solve_lanczos(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, flint_rand_t state, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_solve_block_lanczos(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, slong block_size, flint_rand_t state, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_solve_wiedemann(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_solve_block_wiedemann(gr_ptr x, const gr_lil_mat_t M, gr_srcptr b, slong block_size, flint_rand_t state, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_lil_mat_nullvector_lanczos(gr_ptr x, const gr_lil_mat_t M, flint_rand_t state, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_nullvector_block_lanczos(gr_ptr x, const gr_lil_mat_t M, slong block_size, flint_rand_t state, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_nullvector_wiedemann(gr_ptr x, const gr_lil_mat_t M, flint_rand_t state, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_lil_mat_nullvector_block_wiedemann(gr_ptr x, const gr_lil_mat_t M, slong block_size, flint_rand_t state, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_lil_mat_nullspace(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, const char *algorithm, slong block_size, gr_ctx_t ctx);

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_lil_mat_nullspace_lanczos(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, gr_ctx_t ctx)
{
    return gr_lil_mat_nullspace(X, M, state, max_iters, "lanczos", 1, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_lil_mat_nullspace_wiedemann(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, gr_ctx_t ctx)
{
    return gr_lil_mat_nullspace(X, M, state, max_iters, "wiedemann", 1, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_lil_mat_nullspace_block_lanczos(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, slong block_size, gr_ctx_t ctx)
{
    return gr_lil_mat_nullspace(X, M, state, max_iters, "block lanczos", block_size, ctx);
}

GR_SPARSE_MAT_INLINE WARN_UNUSED_RESULT int gr_lil_mat_nullspace_block_wiedemann(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, slong block_size, gr_ctx_t ctx)
{
    return gr_lil_mat_nullspace(X, M, state, max_iters, "block wiedemann", block_size, ctx);
}


/*
WARN_UNUSED_RESULT int gr_sparse_mat_lu(slong * rank, slong * P, gr_csr_mat_t LU, const gr_csr_mat_t A, int rank_check, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_fflu(slong * res_rank, slong * P, gr_csr_mat_t LU, gr_ptr den, const gr_csr_mat_t A, int rank_check, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_fflu(gr_csr_mat_t X, const gr_csr_mat_t A, const gr_csr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_lu(gr_csr_mat_t X, const gr_csr_mat_t A, const gr_csr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve(gr_csr_mat_t X, const gr_csr_mat_t A, const gr_csr_mat_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_fflu_precomp(gr_csr_mat_t X, const slong * perm, const gr_csr_mat_t A, const gr_csr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_lu_precomp(gr_csr_mat_t X, const slong * perm, const gr_csr_mat_t A, const gr_csr_mat_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_den_fflu(gr_csr_mat_t X, gr_ptr den, const gr_csr_mat_t A, const gr_csr_mat_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_den(gr_csr_mat_t X, gr_ptr den, const gr_csr_mat_t A, const gr_csr_mat_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_solve_field(gr_csr_mat_t X, const gr_csr_mat_t A, const gr_csr_mat_t B, gr_ctx_t ctx);
*/
/*
WARN_UNUSED_RESULT int gr_sparse_mat_det_berkowitz(gr_ptr res, const gr_sparse_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_det_fflu(gr_ptr res, const gr_sparse_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_det_lu(gr_ptr res, const gr_sparse_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_det_cofactor(gr_ptr res, const gr_sparse_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_det_generic_field(gr_ptr res, const gr_sparse_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_det_generic_integral_domain(gr_ptr res, const gr_sparse_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_det_generic(gr_ptr res, const gr_sparse_mat_t A, gr_ctx_t ctx);
*/
/*
WARN_UNUSED_RESULT int gr_sparse_mat_det(gr_ptr res, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_inv(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_adjugate_charpoly(gr_csr_mat_t adj, gr_ptr det, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_adjugate_cofactor(gr_csr_mat_t adj, gr_ptr det, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_adjugate(gr_csr_mat_t adj, gr_ptr det, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_rank_lu(slong * rank, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_rank_fflu(slong * rank, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_rank(slong * rank, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_rref_lu(slong * res_rank, gr_csr_mat_t R, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_rref_fflu(slong * res_rank, gr_csr_mat_t R, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_rref(slong * res_rank, gr_csr_mat_t R, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_rref_den_fflu(slong * res_rank, gr_csr_mat_t R, gr_ptr den, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_rref_den(slong * res_rank, gr_csr_mat_t R, gr_ptr den, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_nullspace(gr_csr_mat_t X, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_ones(gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_pascal(gr_csr_mat_t mat, int triangular, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_stirling(gr_csr_mat_t mat, int kind, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_hilbert(gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_hadamard(gr_csr_mat_t mat, gr_ctx_t ctx);
*/
/* todo: dft, dct */
/*
WARN_UNUSED_RESULT int gr_sparse_mat_transpose(gr_csr_mat_t B, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_tril_classical(gr_csr_mat_t X, const gr_csr_mat_t L, const gr_csr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_tril_recursive(gr_csr_mat_t X, const gr_csr_mat_t L, const gr_csr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_tril(gr_csr_mat_t X, const gr_csr_mat_t L, const gr_csr_mat_t B, int unit, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_triu_classical(gr_csr_mat_t X, const gr_csr_mat_t U, const gr_csr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_triu_recursive(gr_csr_mat_t X, const gr_csr_mat_t U, const gr_csr_mat_t B, int unit, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_nonsingular_solve_triu(gr_csr_mat_t X, const gr_csr_mat_t U, const gr_csr_mat_t B, int unit, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_trace(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly_berkowitz(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_charpoly_berkowitz(gr_poly_t res, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly_danilevsky_inplace(gr_ptr res, gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly_danilevsky(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_charpoly_danilevsky(gr_poly_t res, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly_faddeev(gr_ptr res, gr_csr_mat_t adj, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_charpoly_faddeev(gr_poly_t res, gr_csr_mat_t adj, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly_faddeev_bsgs(gr_ptr res, gr_csr_mat_t adj, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_charpoly_faddeev_bsgs(gr_poly_t res, gr_csr_mat_t adj, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly_from_hessenberg(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_charpoly_from_hessenberg(gr_poly_t cp, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly_gauss(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_charpoly_gauss(gr_poly_t cp, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly_householder(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_charpoly_householder(gr_poly_t cp, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_sparse_mat_charpoly(gr_ptr res, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_charpoly(gr_poly_t res, const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_hessenberg(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_hessenberg_gauss(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_hessenberg_householder(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx);
truth_t gr_sparse_mat_is_hessenberg(const gr_csr_mat_t mat, gr_ctx_t ctx);

int gr_sparse_mat_reduce_row(slong * column, gr_csr_mat_t A, slong * P, slong * L, slong m, gr_ctx_t ctx);
int gr_sparse_mat_apply_row_similarity(gr_csr_mat_t A, slong r, gr_ptr d, gr_ctx_t ctx);
int gr_sparse_mat_minpoly_field(gr_poly_t p, const gr_csr_mat_t X, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_eigenentries(gr_vec_t lambda, gr_vec_t mult, const gr_csr_mat_t mat, int flags, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_eigenentries_other(gr_vec_t lambda, gr_vec_t mult, const gr_csr_mat_t mat, gr_ctx_t mat_ctx, int flags, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_diagonalization_precomp(gr_vec_t D, gr_csr_mat_t L, gr_csr_mat_t R, const gr_csr_mat_t A, const gr_vec_t eigenentries, const gr_vec_t mult, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_diagonalization_generic(gr_vec_t D, gr_csr_mat_t L, gr_csr_mat_t R, const gr_csr_mat_t A, int flags, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_diagonalization(gr_vec_t D, gr_csr_mat_t L, gr_csr_mat_t R, const gr_csr_mat_t A, int flags, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_set_jordan_blocks(gr_csr_mat_t mat, const gr_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_jordan_blocks(gr_vec_t lambda, slong * num_blocks, slong * block_lambda, slong * block_size, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_jordan_transformation(gr_csr_mat_t mat, const gr_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_jordan_form(gr_csr_mat_t J, gr_csr_mat_t P, const gr_csr_mat_t A, gr_ctx_t ctx);

truth_t gr_sparse_mat_is_scalar(const gr_csr_mat_t mat, gr_ctx_t ctx);
truth_t gr_sparse_mat_is_diagonal(const gr_csr_mat_t mat, gr_ctx_t ctx);
truth_t gr_sparse_mat_is_lower_triangular(const gr_csr_mat_t mat, gr_ctx_t ctx);
truth_t gr_sparse_mat_is_upper_triangular(const gr_csr_mat_t mat, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_mul_diag(gr_csr_mat_t C, const gr_csr_mat_t A, const gr_vec_t D, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_diag_mul(gr_csr_mat_t C, const gr_vec_t D, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_exp_jordan(gr_csr_mat_t res, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_exp(gr_csr_mat_t res, const gr_csr_mat_t A, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_mat_log_jordan(gr_csr_mat_t res, const gr_csr_mat_t A, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_mat_log(gr_csr_mat_t res, const gr_csr_mat_t A, gr_ctx_t ctx);
*/

#ifdef __cplusplus
}
#endif

#endif
