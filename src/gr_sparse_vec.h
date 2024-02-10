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

#include <string.h>
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

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_init(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx)
{
    memset(vec, 0, sizeof(gr_sparse_vec_t));
    vec->length = len;
    return GR_SUCCESS;
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_clear(gr_sparse_vec_t vec, gr_ctx_t ctx)
{
    _gr_vec_clear(vec->entries, vec->alloc, ctx);
    flint_free(vec->cols);
    flint_free(vec->entries);
    memset(vec, 0, sizeof(gr_sparse_vec_t));
    return GR_SUCCESS;
}

#define GR_SPARSE_VEC_COL(vec, i) (vec)->cols[i]
#define GR_SPARSE_VEC_ENTRY(vec, i, sz) GR_ENTRY((vec)->entries, i, sz)

GR_SPARSE_VEC_INLINE ulong * gr_sparse_vec_col_ptr(gr_sparse_vec_t vec, slong i, gr_ctx_t ctx) { return vec->cols + i; }
GR_SPARSE_VEC_INLINE const ulong * gr_sparse_vec_col_srcptr(const gr_sparse_vec_t vec, slong i, gr_ctx_t ctx) { return vec->cols + i; }
GR_SPARSE_VEC_INLINE gr_ptr gr_sparse_vec_entry_ptr(gr_sparse_vec_t vec, slong i, gr_ctx_t ctx) { return GR_SPARSE_VEC_ENTRY(vec, i, ctx->sizeof_elem); }
GR_SPARSE_VEC_INLINE gr_srcptr gr_sparse_vec_entry_srcptr(const gr_sparse_vec_t vec, slong i, gr_ctx_t ctx) { return GR_SPARSE_VEC_ENTRY(vec, i, ctx->sizeof_elem); }

void gr_sparse_vec_fit_nnz(gr_sparse_vec_t vec, slong nnz, gr_ctx_t ctx);
void gr_sparse_vec_shrink_to_nnz(gr_sparse_vec_t vec, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE slong gr_sparse_vec_length(const gr_sparse_vec_t vec) { return vec->length; }
void gr_sparse_vec_set_length(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE slong gr_sparse_vec_nnz(const gr_sparse_vec_t vec) { return vec->nnz; }

WARN_UNUSED_RESULT int gr_sparse_vec_set(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_set_from_entries(gr_sparse_vec_t vec, ulong * cols, gr_srcptr entries, slong nnz, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_set_from_entries_sorted_deduped(gr_sparse_vec_t vec, ulong * sorted_deduped_cols, gr_srcptr entries, slong nnz, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_zero(gr_sparse_vec_t vec, gr_ctx_t ctx) { vec->nnz = 0; return GR_SUCCESS; }
WARN_UNUSED_RESULT int gr_sparse_vec_from_dense(gr_sparse_vec_t vec, gr_srcptr src, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_to_dense(gr_ptr vec, gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_slice(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong col_start, slong col_end, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_permute_cols(gr_sparse_vec_t vec, const gr_sparse_vec_t src, slong * p, gr_ctx_t ctx);

truth_t gr_sparse_vec_equal(const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE truth_t gr_sparse_vec_is_zero(const gr_sparse_vec_t vec, gr_ctx_t ctx) { return (vec->nnz == 0 ? T_TRUE : T_FALSE); }

int gr_sparse_vec_write_nz(gr_stream_t out, const gr_sparse_vec_t vec, gr_ctx_t ctx);
int gr_sparse_vec_print_nz(const gr_sparse_vec_t vec, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_randtest(gr_sparse_vec_t vec, double density, slong len, flint_rand_t state, gr_ctx_t ctx);



slong _gr_sparse_vec_count_unique_cols(const ulong *cols0, slong nnz0, const ulong *cols1, slong nnz1);



#define GR_SPV_SWAP_INDS(VEC, I, J, SZ, CTX) \
{                                            \
    slong _temp = (VEC)->cols[I];            \
    (VEC)->cols[I] = (VEC)->cols[J];         \
    (VEC)->cols[J] = _temp;                  \
    gr_swap(GR_ENTRY((VEC)->entries, (I), (SZ)), GR_ENTRY((VEC)->entries, (J), (SZ)), (CTX)); \
}


/* This is used for operations like add which have to riffle through the entries of two vectors */
/* It's an annoying macro that combines with the macros below.  This is so we can handle */
/* Different arguments in different orders, etc */
#define GR_SPV_RFL_TEMPLATE(FUNC_A, FUNC_B, FUNC_AB, DEST_VEC, A_VEC, B_VEC, CTX)                       \
    int status;                                                                                         \
    slong sz, new_nnz, a_ind, b_ind, dest_ind, a_nnz, b_nnz, i;                                         \
    ulong a_col, b_col;                                                                                 \
    if ((DEST_VEC)->length != (A_VEC)->length || (A_VEC)->length != (B_VEC)->length)                    \
        return GR_DOMAIN;                                                                               \
    status = GR_SUCCESS;                                                                                \
    sz = (CTX)->sizeof_elem;                                                                            \
    a_nnz = (A_VEC)->nnz;                                                                               \
    b_nnz = (B_VEC)->nnz;                                                                               \
    new_nnz = _gr_sparse_vec_count_unique_cols((A_VEC)->cols, a_nnz, (B_VEC)->cols, b_nnz);             \
    gr_sparse_vec_fit_nnz((DEST_VEC), new_nnz, (CTX));                                                  \
    /* We go backward through the destination, because it might be an in-place operation on a source */ \
    a_ind = a_nnz-1;                                                                                    \
    b_ind = b_nnz-1;                                                                                    \
    dest_ind = new_nnz-1;                                                                               \
    while (a_ind >= 0 && b_ind >= 0 && status == GR_SUCCESS)                                            \
    {                                                                                                   \
        a_col = (A_VEC)->cols[a_ind];                                                                   \
        b_col = (B_VEC)->cols[b_ind];                                                                   \
        if (a_col > b_col)                                                                              \
        {                                                                                               \
            status |= (FUNC_A);                                                                         \
            (DEST_VEC)->cols[dest_ind] = a_col;                                                         \
            a_ind--;                                                                                    \
        }                                                                                               \
        else if (b_col > a_col)                                                                         \
        {                                                                                               \
            status |= (FUNC_B);                                                                         \
            (DEST_VEC)->cols[dest_ind] = b_col;                                                         \
            b_ind--;                                                                                    \
        }                                                                                               \
        else                                                                                            \
        {                                                                                               \
            status |= (FUNC_AB);                                                                        \
            (DEST_VEC)->cols[dest_ind] = a_col;                                                         \
            a_ind--;                                                                                    \
            b_ind--;                                                                                    \
        }                                                                                               \
        if (T_TRUE != gr_is_zero(GR_ENTRY((DEST_VEC)->entries, dest_ind, sz), (CTX)))                   \
            dest_ind--;                                                                                 \
    }                                                                                                   \
    /* Move the result to the beginning of the dest vec */                                              \
    /* Currently, dest_ind points to one before the start of the legit destination values */            \
    if (dest_ind >= 0 && !status)                                                                       \
    {                                                                                                   \
        new_nnz = (new_nnz-1) - dest_ind;                                                               \
        dest_ind++;                                                                                     \
        for (i = 0; i < new_nnz; i++)                                                                   \
            GR_SPV_SWAP_INDS(DEST_VEC, i, dest_ind + i, sz, CTX);                                       \
    }                                                                                                   \
    (DEST_VEC)->nnz = new_nnz;                                                                          \
    return status;

/* We need some convenience functions */
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_neg_other(gr_ptr res, gr_srcptr src, gr_ctx_t src_ctx, gr_ctx_t ctx) { return (gr_set_other(res, src, src_ctx, ctx) | gr_neg(res, res, ctx)); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_negmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return (gr_mul(res, x, y, ctx) | gr_neg(res, res, ctx)); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_negmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return (gr_mul_si(res, x, y, ctx) | gr_neg(res, res, ctx)); }

#define GR_SPV_RFL_UOP(F, Y, Y_ind) F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, sz), ctx)
#define GR_SPV_RFL_BOP(F, Y, Y_ind, Z, Z_ind) F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, sz), GR_ENTRY(Z->entries, Z_ind, sz), ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_update(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP(gr_set, src, b_ind), GR_SPV_RFL_UOP(gr_set, src, b_ind), res, res, src, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, src1, a_ind), GR_SPV_RFL_UOP(gr_set, src2, b_ind), GR_SPV_RFL_BOP(gr_add, src1, a_ind, src2, b_ind), res, src1, src2, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, src1, a_ind), GR_SPV_RFL_UOP(gr_neg, src2, b_ind), GR_SPV_RFL_BOP(gr_sub, src1, a_ind, src2, b_ind), res, src1, src2, ctx); }

#define GR_SPV_RFL_ZERO gr_zero(GR_ENTRY(res->entries, dest_ind, sz), ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_ZERO, GR_SPV_RFL_ZERO, GR_SPV_RFL_BOP(gr_mul, src1, a_ind, src2, b_ind), res, src1, src2, ctx); }

#define GR_SPV_RFL_UOP_OTHER(F, Y, Y_ind, CTX2) F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, CTX2->sizeof_elem), CTX2, ctx)
#define GR_SPV_RFL_BOP_OTHER(F, Y, Y_ind, Z, Z_ind, CTX2) F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, sz), GR_ENTRY(Z->entries, Z_ind, CTX2->sizeof_elem), CTX2, ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, src1, a_ind), GR_SPV_RFL_UOP_OTHER(gr_set_other, src2, b_ind, ctx2), GR_SPV_RFL_BOP_OTHER(gr_add_other, src1, a_ind, src2, b_ind, ctx2), res, src1, src2, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, src1, a_ind), GR_SPV_RFL_UOP_OTHER(gr_neg_other, src2, a_ind, ctx2), GR_SPV_RFL_BOP_OTHER(gr_sub_other, src1, a_ind, src2, b_ind, ctx2), res, src1, src2, ctx); }

#define GR_SPV_RFL_OTHER_BOP(F, Y, Y_ind, CTX2, Z, Z_ind) F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, CTX2->sizeof_elem), (CTX2), GR_ENTRY(Z->entries, Z_ind, sz), ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_add_sparse_vec(gr_sparse_vec_t res, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP_OTHER(gr_set_other, src1, a_ind, ctx1), GR_SPV_RFL_UOP(gr_set, src2, b_ind), GR_SPV_RFL_OTHER_BOP(gr_other_add, src1, a_ind, ctx1, src2, b_ind), res, src1, src2, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_sub_sparse_vec(gr_sparse_vec_t res, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP_OTHER(gr_set_other, src1, a_ind, ctx1), GR_SPV_RFL_UOP(gr_neg, src2, b_ind), GR_SPV_RFL_OTHER_BOP(gr_other_sub, src1, a_ind, ctx1, src2, b_ind), res, src1, src2, ctx); }

#define GR_SPV_RFL_UOP_SCALAR(F, Y, Y_ind) F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, sz), c, ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP_SCALAR(gr_mul, src, b_ind), GR_SPV_RFL_UOP_SCALAR(gr_addmul, src, b_ind), res, res, src, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP_SCALAR(gr_negmul, src, b_ind), GR_SPV_RFL_UOP_SCALAR(gr_submul, src, b_ind), res, res, src, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP_SCALAR(gr_mul_si, src, b_ind), GR_SPV_RFL_UOP_SCALAR(gr_addmul_si, src, b_ind), res, res, src, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP_SCALAR(gr_negmul_si, src, b_ind), GR_SPV_RFL_UOP_SCALAR(gr_submul_si, src, b_ind), res, res, src, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_neg(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx) { return (gr_sparse_vec_set(res, src, ctx) | _gr_vec_neg(res->entries, res->entries, res->nnz, ctx)); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sum(gr_ptr res, const gr_sparse_vec_t vec, gr_ctx_t ctx) { return _gr_vec_sum(res, vec->entries, vec->nnz, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_nz_product(gr_ptr res, const gr_sparse_vec_t vec, gr_ctx_t ctx) { return _gr_vec_product(res, vec->entries, vec->nnz, ctx); }


/* These are versions that incorporate a sparse into a dense */
#define GR_SPV_INTO_DENSE_TEMPLATE(FUNC, SVEC, CTX) \
    slong i;                                        \
    slong sz = (CTX)->sizeof_elem;                  \
    int status = GR_SUCCESS;                        \
    slong nnz = (SVEC)->nnz;                        \
    for (i = 0; i < nnz; i++)                       \
    {                                               \
        status |= (FUNC);                           \
    }                                               \
    return status;

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_update_to_dense(gr_ptr dres, const gr_sparse_vec_t src, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_set(GR_ENTRY(dres, src->cols[i], sz), GR_ENTRY(src->entries, i, sz), ctx), src, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_to_dense(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_add(GR_ENTRY(dres, svec2->cols[i], sz), GR_ENTRY(dvec1, svec2->cols[i], sz), GR_ENTRY(svec2->entries, i, sz), ctx), svec2, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_to_dense(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_sub(GR_ENTRY(dres, svec2->cols[i], sz), GR_ENTRY(dvec1, svec2->cols[i], sz), GR_ENTRY(svec2->entries, i, sz), ctx), svec2, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_nz_mul_to_dense(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_mul(GR_ENTRY(dres, svec2->cols[i], sz), GR_ENTRY(dvec1, svec2->cols[i], sz), GR_ENTRY(svec2->entries, i, sz), ctx), svec2, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_nz_div_to_dense(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_div(GR_ENTRY(dres, svec2->cols[i], sz), GR_ENTRY(dvec1, svec2->cols[i], sz), GR_ENTRY(svec2->entries, i, sz), ctx), svec2, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar_to_dense(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_addmul(GR_ENTRY(dres, svec->cols[i], sz), GR_ENTRY(svec->entries, i, sz), c, ctx), svec, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar_to_dense(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_submul(GR_ENTRY(dres, svec->cols[i], sz), GR_ENTRY(svec->entries, i, sz), c, ctx), svec, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar_si_to_dense(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_addmul_si(GR_ENTRY(dres, svec->cols[i], sz), GR_ENTRY(svec->entries, i, sz), c, ctx), svec, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar_si_to_dense(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_submul_si(GR_ENTRY(dres, svec->cols[i], sz), GR_ENTRY(svec->entries, i, sz), c, ctx), svec, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar_fmpz_to_dense(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_addmul_fmpz(GR_ENTRY(dres, svec->cols[i], sz), GR_ENTRY(svec->entries, i, sz), c, ctx), svec, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar_fmpz_to_dense(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx) { GR_SPV_INTO_DENSE_TEMPLATE(gr_submul_fmpz(GR_ENTRY(dres, svec->cols[i], sz), GR_ENTRY(svec->entries, i, sz), c, ctx), svec, ctx) }

/***** TODO: 

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_OP(ctx, VEC_DOT)(res, initial, subtract, vec1, vec2, len, ctx); }

*/




#define GR_SPARSE_VEC_DENSE_VEC_OP(dense_vec_op, res, src, c, ctx)     \
    if(res->length != src->length)                                     \
    {                                                                  \
        return GR_DOMAIN;                                              \
    }                                                                  \
    if(res != src)                                                     \
    {                                                                  \
        gr_sparse_vec_fit_nnz(res, src->nnz, ctx);                     \
        res->nnz = src->nnz;                                           \
        memcpy(res->cols, src->cols, src->nnz*sizeof(slong));          \
    }                                                                  \
    return dense_vec_op(res->entries, src->entries, src->nnz, c, ctx); \

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar, res, src, c, ctx) } 
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_si, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_ui, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpz, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpq, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_2exp_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_2exp_si, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar, res, src, c, ctx) } 
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_si, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_ui, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_fmpz, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_fmpq, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar, res, src, c, ctx) } 
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_si, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_ui, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpz, res, src, c, ctx) }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpq, res, src, c, ctx) }



#ifdef __cplusplus
}
#endif

#endif
