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
    ulong *cols;
    slong length;
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
    flint_free(vec->cols);
    flint_free(vec->entries);
    memset(vec, 0, sizeof(gr_sparse_vec_t));
    return GR_SUCCESS;
}

#define GR_SPARSE_VEC_ENTRY(vec, i, sz) GR_ENTRY((vec)->entries, i, sz)

GR_SPARSE_VEC_INLINE ulong *
gr_sparse_vec_col_ptr(gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
{
    return vec->cols + i;
}

GR_SPARSE_VEC_INLINE const ulong *
gr_sparse_vec_col_srcptr(const gr_sparse_vec_t vec, slong i, gr_ctx_t ctx)
{
    return vec->cols + i;
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
slong gr_sparse_vec_length(const gr_sparse_vec_t vec) { return vec->length; }
void gr_sparse_vec_set_length(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_set(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx);

slong _gr_sparse_vec_count_unique_cols(const ulong *cols0, slong nnz0, const ulong *cols1, slong nnz1);

/* This is used for operations like add */
#define GR_SPARSE_VEC_RIFFLE_TEMPLATE(LIT_FUNC_A, LIT_FUNC_B, LIT_FUNC_AB, DEST_VEC, A_COLS, A_NNZ, B_COLS, B_NNZ, CTX, CTX2) \
{
    int status = GR_SUCCESS;
    slong sz = (CTX)->sizeof_elem;
    slong sz2 = (CTX2)->sizeof_elem;
    slong new_nnz = _gr_sparse_vec_count_unique_cols((A_COLS), (A_NNZ), (B_COLS), (B_NNZ))
    gr_sparse_vec_fit_nnz(DEST_VEC, new_nnz, (CTX));
    /* We go backward through the destination, because it might be an in-place operation on a source */
    slong a_ind = (A_NNZ)-1;
    slong b_ind = (B_NNZ)-1;
    slong dst_ind = new_nnz-1;
    while (a_ind >= 0 && b_ind >= 0 && status == GR_SUCCESS)
    {
        slong a_col = (A_COLS)[a_ind];
        slong b_col = (B_COLS)[b_ind];
        if (a_col > b_col)
        {
            status |= LIT_FUNC_A;
            (DEST_VEC)->cols[dst_ind] = a_col;
            a_ind--;
        }
        else if (b_col > a_col)
        {
            status |= LIT_FUNC_B;
            (DEST_VEC)->cols[dst_ind] = b_col;
            b_ind--;
        }
        else
        {
            status |= LIT_FUNC_AB;
            (DEST_VEC)->cols[dest_ind] = a_col;
            a_ind--;
            b_ind--;
        }
        dst_ind--;
    }
    (DEST_VEC)->nnz = new_nnz;
    return status;
}

#define GR_SVEC_RIFFLE_COPY(Y, Y_ind) gr_set(GR_ENTRY(res->entries, dst_ptr, sz), GR_ENTRY(Y->entries, Y_ptr, sz), ctx)
#define GR_SVEC_RIFFLE_NEG(Y, Y_ind) gr_neg(GR_ENTRY(res->entries, dst_ptr, sz), GR_ENTRY(Y->entries, Y_ptr, sz), ctx)
#define GR_SVEC_RIFFLE_ADD() gr_add(GR_ENTRY(res->entries, dst_ptr, sz), GR_ENTRY(src1->entries, a_ptr, sz), GR_ENTRY(src2->entries, b_ptr, sz), ctx)
#define SUBAB gr_sub(GR_ENTRY(res->entries, dst_ptr, sz), GR_ENTRY(src1->entries, a_ptr, sz), GR_ENTRY(src2->entries, b_ptr, sz), ctx)
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_update(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx) { GR_SPARSE_VEC_RIFFLE_TEMPLATE(COPYRES, COPYSRC, COPYSRC, res, res->inds, res->nnz, src->inds, src->nnz, ctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { GR_SPARSE_VEC_RIFFLE_TEMPLATE(COPYSRC1, COPYSRC2, ADDAB, vec, src1->inds, src1->nnz, src2->inds, src2->nnz, ctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { GR_SPARSE_VEC_RIFFLE_TEMPLATE(COPYSRC1, NEGSRC2, SUBAB, vec, src1->inds, src1->nnz, src2->inds, src2->nnz, ctx, ctx); }

#define COPYOTHERSRC2 gr_set_other(GR_ENTRY(res->entries, dst_ptr, sz), GR_ENTRY(src2->entries, b_ptr, sz), ctx2, ctx)
#define ADDOTHERAB gr_add_other(GR_ENTRY(res->entries, dst_ptr, sz), GR_ENTRY(src1->entries, a_ptr, sz), GR_ENTRY(src2->entries, b_ptr, sz), ctx2, ctx)
#define SUBOTHERAB gr_sub_other(GR_ENTRY(res->entries, dst_ptr, sz), GR_ENTRY(src1->entries, a_ptr, sz), GR_ENTRY(src2->entries, b_ptr, sz), ctx2, ctx)
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, slong len, gr_ctx_t ctx) { GR_SPARSE_VEC_RIFFLE_TEMPLATE(COPYSRC1, COPYOTHERSRC2, ADDAOTHERB, vec, src1->inds, src1->nnz, src2->inds, src2->nnz, ctx, ctx2); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, slong len, gr_ctx_t ctx) { GR_SPARSE_VEC_RIFFLE_TEMPLATE(COPYSRC1, COPYOTHERSRC2, SUBAOTHERB, vec, src1->inds, src1->nnz, src2->inds, src2->nnz, ctx, ctx2); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_other_add(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx2, const gr_sparse_vec_t vec3, slong len, gr_ctx_t ctx) { return GR_OTHER_OP_VEC(ctx, OTHER_ADD_VEC)(vec1, vec2, ctx2, vec3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_other_sub(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx2, const gr_sparse_vec_t vec3, slong len, gr_ctx_t ctx) { return GR_OTHER_OP_VEC(ctx, OTHER_SUB_VEC)(vec1, vec2, ctx2, vec3, len, ctx); }



#define GR_SPARSE_VEC_DENSE_VEC_OP(g, res, src, c, ctx) \
{  \
    int status = GR_SUCCESS; \
    status = gr_vec_set(res, src, ctx); \
    if(status == GR_SUCCESS) { \
      status = g(res->entries, src->entries, len, c, ctx); \
    } \
    return status; \
}

int gr_sparse_vec_mul_scalar(gr_vec_t res, gr_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar, blah, src, c, ctx) } 
int gr_sparse_vec_mul_scalar_si(gr_vec_t res, gr_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_si, res, src, c, ctx) }





//I think the below might be the same as "add"
WARN_UNUSED_RESULT int gr_sparse_vec_extend(gr_sparse_vec_t vec, const gr_sparse_vec_t src, gr_ctx_t ctx);

int gr_sparse_vec_write(gr_stream_t out, const gr_sparse_vec_t vec, gr_ctx_t ctx);
int gr_sparse_vec_print(const gr_sparse_vec_t vec, gr_ctx_t ctx);

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_zero(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_CONSTANT_OP(ctx, VEC_ZERO)(vec, len, ctx); }
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


GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_MUL)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_DIV)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_DIVEXACT)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_POW)(res, src1, src2, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_MUL_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_DIV_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_DIVEXACT_SCALAR)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP(ctx, VEC_POW_SCALAR)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_MUL_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_DIV_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_DIVEXACT_SCALAR_SI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_POW_SCALAR_SI)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_MUL_SCALAR_UI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_DIV_SCALAR_UI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_DIVEXACT_SCALAR_UI)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_ui(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, ulong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_UI(ctx, VEC_POW_SCALAR_UI)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_MUL_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_DIV_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_DIVEXACT_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_fmpz(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpz_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPZ(ctx, VEC_POW_SCALAR_FMPZ)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_MUL_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_DIV_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_DIVEXACT_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_fmpq(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, const fmpq_t c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_FMPQ(ctx, VEC_POW_SCALAR_FMPQ)(vec1, vec2, len, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_mul_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_MUL_VEC)(vec1, c, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_div_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_DIV_VEC)(vec1, c, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_divexact_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_DIVEXACT_VEC)(vec1, c, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_pow_vec(gr_sparse_vec_t vec1, gr_srcptr c, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_VEC_OP(ctx, SCALAR_POW_VEC)(vec1, c, vec2, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_MUL_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_DIV_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_DIVEXACT_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, const gr_sparse_vec_t vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_OTHER(ctx, VEC_POW_OTHER)(vec1, vec2, vec3, ctx3, len, ctx); }

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
