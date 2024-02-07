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
#define GR_SPARSE_VEC_INLINE static __inline__
#endif

#include "gr.h"
#include "gr_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    slong alloc;
    gr_ptr entries;
    ulong *inds;
    slong len;
    slong nnz;
}
gr_sparse_vec_struct;

typedef gr_sparse_vec_struct gr_sparse_vec_t[1];

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_init(gr_sparse_vec_t vec, gr_ctx_t ctx) {
    memset(vec, 0, sizeof(gr_sparse_vec_t));
    return GR_SUCCESS;
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_clear(gr_sparse_vec_t vec, gr_ctx_t ctx) {
    _gr_vec_clear(vec->entries, vec->alloc, ctx);
    flint_free(vec->inds);
    flint_free(vec->entries);
    memset(vec, 0, sizeof(gr_sparse_vec_t));
    return GR_SUCCESS;
}

#define GR_SPARSE_VEC_ENTRY(vec, i, sz) GR_ENTRY((vec)->entries, i, sz)

GR_SPARSE_VEC_INLINE ulong *
gr_sparse_vec_ind_ptr(gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
{
    return vec->inds + i;
}

GR_SPARSE_VEC_INLINE const ulong *
gr_sparse_vec_ind_srcptr(const gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
{
    return vec->inds + i;
}

GR_SPARSE_VEC_INLINE gr_ptr
gr_sparse_vec_entry_ptr(gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
{
    return GR_SPARSE_VEC_ENTRY(vec, i, ctx->sizeof_elem);
}

GR_SPARSE_VEC_INLINE gr_srcptr
gr_sparse_vec_entry_srcptr(const gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
{
    return GR_SPARSE_VEC_ENTRY(vec, i, ctx->sizeof_elem);
}

GR_SPARSE_VEC_INLINE slong gr_sparse_vec_nnz(const gr_sparse_vec_t vec, gr_ctx_t ctx)
{
    return vec->nnz;
}

void gr_sparse_vec_fit_nnz(gr_sparse_vec_t vec, slong nnz, gr_ctx_t ctx);
void gr_sparse_vec_set_length(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_set(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_update(gr_sparse_vec_t vec, const gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_update_entries(gr_sparse_vec_t vec, gr_srcptr entries, const ulong* inds, slong nnz, gr_ctx_t ctx);
//I think the below might be the same as "add"
WARN_UNUSED_RESULT int gr_sparse_vec_extend(gr_sparse_vec_t vec, const gr_sparse_vec_t src, gr_ctx_t ctx);

int gr_sparse_vec_write(gr_stream_t out, const gr_sparse_vec_t vec, gr_ctx_t ctx);
int gr_sparse_vec_print(const gr_sparse_vec_t vec, gr_ctx_t ctx);

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_zero(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_CONSTANT_OP(ctx, VEC_ZERO)(vec, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_set(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP(ctx, VEC_SET)(res, src, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_neg(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP(ctx, VEC_NEG)(res, src, len, ctx); }

WARN_UNUSED_RESULT int gr_sparse_vec_from_entries(gr_sparse_vec_t vec, slong * inds, gr_srcptr vals, slong nnz);

WARN_UNUSED_RESULT int gr_sparse_vec_from_dense(gr_sparse_vec_t vec, gr_srcptr src, slong len);

WARN_UNUSED_RESULT int gr_sparse_vec_to_dense(gr_ptr vec, gr_sparse_vec_t src, slong len);
WARN_UNUSED_RESULT int gr_sparse_vec_split(gr_sparse_vec_t res1, gr_sparse_vec_t res2, const gr_sparse_vec_t vec, slong ind);

/* Vector permutation */
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_permute_inds(gr_sparse_vec_t vec, slong *P) 
{
    slong i;
    for (i = 0; i < vec->nnz; ++i) vec->inds[i] = P[vec->inds[i]];
    // TODO: re-sort vector
}

/*
typedef int ((*gr_method_vec_normalise_op)(slong *, gr_srcptr, slong, gr_ctx_ptr));
typedef slong ((*gr_method_vec_normalise_weak_op)(gr_srcptr, slong, gr_ctx_ptr));
#define GR_SPARSE_VEC_NORMALISE_OP(ctx, NAME) (((gr_method_vec_normalise_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_VEC_NORMALISE_WEAK_OP(ctx, NAME) (((gr_method_vec_normalise_weak_op *) ctx->methods)[GR_METHOD_ ## NAME])

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_normalise(slong * res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_NORMALISE_OP(ctx, VEC_NORMALISE)(res, vec, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT slong gr_sparse_vec_normalise_weak(const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_NORMALISE_WEAK_OP(ctx, VEC_NORMALISE_WEAK)(vec, len, ctx); }
*/

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_ADD)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_SUB)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_MUL)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_DIV)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_DIVEXACT)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_POW)(res, src1, src2, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_ADD_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_SUB_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_MUL_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_DIV_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_DIVEXACT_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_POW_SCALAR)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_ADD_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_SUB_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_MUL_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_DIV_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_DIVEXACT_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_POW_SCALAR_SI)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_ADD_SCALAR_UI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_SUB_SCALAR_UI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_MUL_SCALAR_UI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_DIV_SCALAR_UI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_DIVEXACT_SCALAR_UI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_POW_SCALAR_UI)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_ADD_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_SUB_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_MUL_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_DIV_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_DIVEXACT_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_POW_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_ADD_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_SUB_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_MUL_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_DIV_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_DIVEXACT_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_POW_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_add_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_ADD_VEC)(vec1, c, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_sub_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_SUB_VEC)(vec1, c, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_mul_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_MUL_VEC)(vec1, c, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_div_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_DIV_VEC)(vec1, c, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_divexact_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_DIVEXACT_VEC)(vec1, c, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_pow_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_POW_VEC)(vec1, c, vec2, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_ADD_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_SUB_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_MUL_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_DIV_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_DIVEXACT_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_POW_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_add_vec(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx2, const gr_sparse_vec_t vec3, slong len, gr_ctx_t ctx) { return GR_OTHER_OP_VEC(ctx, OTHER_ADD_VEC)(vec1, vec2, ctx2, vec3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_sub_vec(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx2, const gr_sparse_vec_t vec3, slong len, gr_ctx_t ctx) { return GR_OTHER_OP_VEC(ctx, OTHER_SUB_VEC)(vec1, vec2, ctx2, vec3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_mul_vec(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx2, const gr_sparse_vec_t vec3, slong len, gr_ctx_t ctx) { return GR_OTHER_OP_VEC(ctx, OTHER_MUL_VEC)(vec1, vec2, ctx2, vec3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_div_vec(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx2, const gr_sparse_vec_t vec3, slong len, gr_ctx_t ctx) { return GR_OTHER_OP_VEC(ctx, OTHER_DIV_VEC)(vec1, vec2, ctx2, vec3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_divexact_vec(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx2, const gr_sparse_vec_t vec3, slong len, gr_ctx_t ctx) { return GR_OTHER_OP_VEC(ctx, OTHER_DIVEXACT_VEC)(vec1, vec2, ctx2, vec3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_pow_vec(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx2, const gr_sparse_vec_t vec3, slong len, gr_ctx_t ctx) { return GR_OTHER_OP_VEC(ctx, OTHER_POW_VEC)(vec1, vec2, ctx2, vec3, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_ADD_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_SUB_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_MUL_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_DIV_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_DIVEXACT_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_POW_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_add_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_ADD_VEC)(vec1, c, cctx, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_sub_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_SUB_VEC)(vec1, c, cctx, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_mul_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_MUL_VEC)(vec1, c, cctx, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_div_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_DIV_VEC)(vec1, c, cctx, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_divexact_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_DIVEXACT_VEC)(vec1, c, cctx, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_pow_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_POW_VEC)(vec1, c, cctx, vec2, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_2exp_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_MUL_SCALAR_2EXP_SI)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_ADDMUL_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_SUBMUL_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_ADDMUL_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_SUBMUL_SCALAR_SI)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE truth_t gr_sparse_vec_equal(const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_PREDICATE(ctx, VEC_EQUAL)(vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE truth_t gr_sparse_vec_is_zero(const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_PREDICATE(ctx, VEC_IS_ZERO)(vec, len, ctx); }

typedef int ((*gr_method_vec_reduce_op)(gr_sparse_vec_t, gr_srcptr, slong, gr_ctx_ptr));

typedef int ((*gr_method_vec_dot_op)(gr_sparse_vec_t, gr_srcptr, int, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_dot_si_op)(gr_sparse_vec_t, gr_srcptr, int, gr_srcptr, const slong *, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_dot_ui_op)(gr_sparse_vec_t, gr_srcptr, int, gr_srcptr, const ulong *, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_dot_fmpz_op)(gr_sparse_vec_t, gr_srcptr, int, gr_srcptr, const fmpz *, slong, gr_ctx_ptr));

#define GR_SPARSE_VEC_REDUCE_OP(ctx, NAME) (((gr_method_vec_reduce_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_VEC_DOT_OP(ctx, NAME) (((gr_method_vec_dot_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_VEC_DOT_SI_OP(ctx, NAME) (((gr_method_vec_dot_si_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_VEC_DOT_UI_OP(ctx, NAME) (((gr_method_vec_dot_ui_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SPARSE_VEC_DOT_FMPZ_OP(ctx, NAME) (((gr_method_vec_dot_fmpz_op *) ctx->methods)[GR_METHOD_ ## NAME])

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sum(gr_ptr res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_REDUCE_OP(ctx, VEC_SUM)(res, vec, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_product(gr_ptr res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_REDUCE_OP(ctx, VEC_PRODUCT)(res, vec, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_OP(ctx, VEC_DOT)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_OP(ctx, VEC_DOT_REV)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const slong * vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_SI_OP(ctx, VEC_DOT_SI)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const ulong * vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_UI_OP(ctx, VEC_DOT_UI)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const fmpz * vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_FMPZ_OP(ctx, VEC_DOT_FMPZ)(res, initial, subtract, vec1, vec2, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_reciprocals(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, VEC_RECIPROCALS)(res, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP(ctx, VEC_SET_POWERS)(res, x, len, ctx); }

/* todo: could allow overloading this as well */
/* todo: worth warning about unused result? */
WARN_UNUSED_RESULT int gr_sparse_vec_randtest(gr_ptr res, flint_rand_t state, slong len, gr_ctx_t ctx);

/*
WARN_UNUSED_RESULT int gr_sparse_vec_sum_bsplit_parallel(gr_ptr res, const gr_sparse_vec_t vec, slong len, slong basecase_cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_sum_bsplit(gr_ptr res, const gr_sparse_vec_t vec, slong len, slong basecase_cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_sum_parallel(gr_ptr res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_sum_serial(gr_ptr res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_sum_generic(gr_ptr res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_product_bsplit_parallel(gr_ptr res, const gr_sparse_vec_t vec, slong len, slong basecase_cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_product_bsplit(gr_ptr res, const gr_sparse_vec_t vec, slong len, slong basecase_cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_product_parallel(gr_ptr res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_product_serial(gr_ptr res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_product_generic(gr_ptr res, const gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_step(gr_sparse_vec_t vec, gr_srcptr start, gr_srcptr step, slong len, gr_ctx_t ctx);
*/

#ifdef __cplusplus
}
#endif

#endif
