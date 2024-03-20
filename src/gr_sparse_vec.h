/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef GR_SPARSE_VEC_H
#define GR_SPARSE_VEC_H

#ifdef GR_SPARSE_VEC_INLINES_C
#define GR_SPARSE_VEC_INLINE
#else
#define GR_SPARSE_VEC_INLINE static inline
#endif

#include <string.h>
#include "gr.h"
#include "gr_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

/**
 * Types and basic access
**/
typedef struct
{
    slong length;
    slong nnz;
    slong alloc;
    ulong *inds;
    gr_ptr nzs;
}
gr_sparse_vec_struct;

typedef gr_sparse_vec_struct gr_sparse_vec_t[1];

GR_SPARSE_VEC_INLINE void
gr_sparse_vec_init(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx)
{
    memset(vec, 0, sizeof(gr_sparse_vec_t));
    vec->length = len;
}

GR_SPARSE_VEC_INLINE void
gr_sparse_vec_clear(gr_sparse_vec_t vec, gr_ctx_t ctx)
{
    _gr_vec_clear(vec->nzs, vec->alloc, ctx);
    flint_free(vec->inds);
    flint_free(vec->nzs);
    memset(vec, 0, sizeof(gr_sparse_vec_t));
}

GR_SPARSE_VEC_INLINE slong 
gr_sparse_vec_length(const gr_sparse_vec_t vec) 
{ return vec->length; }

void gr_sparse_vec_set_length(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);

GR_SPARSE_VEC_INLINE slong 
gr_sparse_vec_nnz(const gr_sparse_vec_t vec)
{ return vec->nnz; }

void gr_sparse_vec_fit_nnz(gr_sparse_vec_t vec, slong nnz, gr_ctx_t ctx);
void gr_sparse_vec_shrink_to_nnz(gr_sparse_vec_t vec, gr_ctx_t ctx);

WARN_UNUSED_RESULT int
gr_sparse_vec_from_entries(gr_sparse_vec_t vec, ulong * inds, gr_srcptr entries, slong nnz, int is_canonical, gr_ctx_t ctx);

WARN_UNUSED_RESULT int 
gr_sparse_vec_randtest(gr_sparse_vec_t vec, slong nnz, int replacement, flint_rand_t state, gr_ctx_t ctx);

WARN_UNUSED_RESULT int 
gr_sparse_vec_randtest_prob(gr_sparse_vec_t vec, double prob, flint_rand_t state, gr_ctx_t ctx);

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_is_valid(gr_sparse_vec_t vec, gr_ctx_t ctx)
{
    slong i, sz = ctx->sizeof_elem;

    // Check that parameters are valid
    if (vec->nnz > vec->alloc || vec->alloc > vec->length)
        return 0;

    // Check that entries are valid
    for (i = 0; i < vec->nnz; ++i)
    {
        if (vec->inds[i] >= vec->length || (i > 0 && vec->inds[i] <= vec->inds[i-1]))
            return 0;
        if (gr_is_zero(GR_ENTRY(vec->nzs, i, sz), ctx) == T_TRUE)
            return 0;
    }
    return 1;
}

/**
 * Getting, setting, and conversion
**/

#define GR_SPARSE_VEC_IND(vec, nz_idx) (vec)->inds[nz_idx]
#define GR_SPARSE_VEC_ENTRY(vec, nz_idx, sz) GR_ENTRY((vec)->nzs, nz_idx, sz)

GR_SPARSE_VEC_INLINE ulong * 
gr_sparse_vec_ind_ptr(gr_sparse_vec_t vec, slong nz_idx, gr_ctx_t ctx)
{
    if (nz_idx < 0 || nz_idx >= vec->nnz)
        return NULL;
    return vec->inds + nz_idx;
}

GR_SPARSE_VEC_INLINE const ulong * 
gr_sparse_vec_ind_srcptr(const gr_sparse_vec_t vec, slong nz_idx, gr_ctx_t ctx)
{
    if (nz_idx < 0 || nz_idx >= vec->nnz)
        return NULL;
    return vec->inds + nz_idx;
}

GR_SPARSE_VEC_INLINE gr_ptr 
gr_sparse_vec_entry_ptr(gr_sparse_vec_t vec, slong nz_idx, gr_ctx_t ctx)
{
    if (nz_idx < 0 || nz_idx >= vec->nnz)
        return NULL;
    return GR_SPARSE_VEC_ENTRY(vec, nz_idx, ctx->sizeof_elem);
}

GR_SPARSE_VEC_INLINE gr_srcptr 
gr_sparse_vec_entry_srcptr(const gr_sparse_vec_t vec, slong nz_idx, gr_ctx_t ctx)
{
    if (nz_idx < 0 || nz_idx >= vec->nnz)
        return NULL;
    return GR_SPARSE_VEC_ENTRY(vec, nz_idx, ctx->sizeof_elem);
}

WARN_UNUSED_RESULT gr_ptr gr_sparse_vec_find_entry(gr_sparse_vec_t vec, slong ind, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_get_entry(gr_ptr dst, gr_sparse_vec_t vec, slong ind, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_set_entry(gr_sparse_vec_t vec, slong ind, gr_srcptr entry, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_set(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_slice(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong ind_start, slong ind_end, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_set_vec(gr_sparse_vec_t dst, gr_srcptr src, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_set_sparse_vec(gr_ptr dst, gr_sparse_vec_t src, gr_ctx_t ctx);

GR_SPARSE_VEC_INLINE void gr_sparse_vec_swap(gr_sparse_vec_t vec1, gr_sparse_vec_t vec2, gr_ctx_t ctx)
{
    FLINT_SWAP(slong, vec1->alloc, vec2->alloc);
    FLINT_SWAP(gr_ptr, vec1->nzs, vec2->nzs);
    FLINT_SWAP(ulong *, vec1->inds, vec2->inds);
    FLINT_SWAP(slong, vec1->length, vec2->length);
    FLINT_SWAP(slong, vec1->nnz, vec2->nnz);
}

GR_SPARSE_VEC_INLINE 
void gr_sparse_vec_zero(gr_sparse_vec_t vec, gr_ctx_t ctx)
{ vec->nnz = 0; }

GR_SPARSE_VEC_INLINE 
int gr_sparse_vec_one(gr_sparse_vec_t vec, slong ind, gr_ctx_t ctx)
{
    if (ind < 0)
        return GR_DOMAIN;

    gr_sparse_vec_fit_nnz(vec, 1, ctx);
    vec->inds[0] = ind;
    vec->nnz = 1;
    return gr_one(vec->nzs, ctx);
}

WARN_UNUSED_RESULT int gr_sparse_vec_permute_inds(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong * p, gr_ctx_t ctx);

/**
 * Comparison
**/

truth_t gr_sparse_vec_equal(const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE truth_t gr_sparse_vec_is_zero(const gr_sparse_vec_t vec, gr_ctx_t ctx) { return _gr_vec_is_zero(vec->nzs, vec->nnz, ctx); }

/**
 * Output
**/

int gr_sparse_vec_write_nz(gr_stream_t out, const gr_sparse_vec_t vec, gr_ctx_t ctx);
int gr_sparse_vec_print_nz(const gr_sparse_vec_t vec, gr_ctx_t ctx);

/**
 * Arithmetic
**/

// Internal function to count the union of indicies in inds0 and inds1
slong _gr_sparse_vec_count_unique_inds(const ulong *inds0, slong nnz0, const ulong *inds1, slong nnz1);

WARN_UNUSED_RESULT int gr_sparse_vec_neg(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_update(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_add(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_sub(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx) ;
WARN_UNUSED_RESULT int gr_sparse_vec_mul(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_add_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_sub_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_mul_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_other_add_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_other_sub_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_other_mul_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx);

/**
 * Arithmetic into dense vectors
**/

WARN_UNUSED_RESULT int gr_vec_update_sparse_vec_nz(gr_ptr dres, const gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_add_sparse_vec(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_sub_sparse_vec(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_mul_sparse_vec_nz(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_div_sparse_vec_nz(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_vec_addmul_sparse_vec_scalar(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_submul_sparse_vec_scalar(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_addmul_sparse_vec_scalar_si(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_submul_sparse_vec_scalar_si(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_addmul_sparse_vec_scalar_fmpz(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_submul_sparse_vec_scalar_fmpz(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx);

/**
 * Scalar multiplication and division
**/

WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_2exp_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx);

/**
 * Sum and product
**/

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_sum(gr_ptr dst, const gr_sparse_vec_t src, gr_ctx_t ctx)
{ return _gr_vec_sum(dst, src->nzs, src->nnz, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_nz_product(gr_ptr dst, const gr_sparse_vec_t src, gr_ctx_t ctx)
{ return _gr_vec_product(dst, src->nzs, src->nnz, ctx); }

/**
 * Dot product
**/

WARN_UNUSED_RESULT int gr_sparse_vec_dot(gr_ptr dst, gr_srcptr initial, int subtract, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_dot_vec(gr_ptr dst, gr_srcptr initial, int subtract, const gr_sparse_vec_t src1, gr_srcptr src2, gr_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
