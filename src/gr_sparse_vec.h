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

WARN_UNUSED_RESULT int gr_sparse_vec_set_from_entries(gr_sparse_vec_t vec, slong * cols, gr_srcptr entries, slong nnz);
WARN_UNUSED_RESULT int gr_sparse_vec_set_from_entries_sorted_deduped(gr_sparse_vec_t vec, slong * sorted_cols, gr_srcptr entries, slong nnz);
WARN_UNUSED_RESULT int gr_sparse_vec_from_dense(gr_sparse_vec_t vec, gr_srcptr src, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_to_dense(gr_ptr vec, gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_slice(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong col_start, slong col_end, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_permute_cols(gr_sparse_vec_t vec, const gr_sparse_vec_t src, slong * p, gr_ctx_t ctx);

truth_t gr_sparse_vec_equal(const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE truth_t gr_sparse_vec_is_zero(const gr_sparse_vec_t vec, gr_ctx_t ctx) { return (vec->nnz == 0 ? T_TRUE ? T_FALSE); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_zero(gr_sparse_vec_t vec, gr_ctx_t ctx) { vec->nnz = 0; return GR_SUCCESS; }

int gr_sparse_vec_write_nz(gr_stream_t out, const gr_sparse_vec_t vec, gr_ctx_t ctx);
int gr_sparse_vec_print_nz(const gr_sparse_vec_t vec, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_randtest(gr_sparse_vec_t vec, double density, slong len, flint_rand_t state, gr_ctx_t ctx);



slong _gr_sparse_vec_count_unique_cols(const ulong *cols0, slong nnz0, const ulong *cols1, slong nnz1);



#define GR_SPV_SWAP_INDS(VEC, I, J, CTX) \
{                                        \
    slong _temp = (VEC)->cols[I];        \
    (VEC)->cols[I] = (VEC)->cols[J];     \
    (VEC)->cols[J] = _temp;              \
    gr_swap(GR_ENTRY((VEC)->entries, (I), ctx), GR_ENTRY((VEC)->entries, (J), ctx), ctx); \
}


/* This is used for operations like add which have to riffle through the entries of two vectors */
#define GR_SPV_RFL_TEMPLATE(FUNC_A, FUNC_B, FUNC_AB, DEST_VEC, A_VEC, B_VEC, CTX)                       \
    int status;                                                                                         \
    slong sz, new_nnz, a_ind, b_ind, dest_ind, a_nnz, b_nnz, a_col, b_col, i;                           \
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
        if (T_FALSE == gr_is_zero(GR_ENTRY((DEST_VEC)->entries, dest_ind, sz)))                         \
            dest_ind--;                                                                                 \
    }                                                                                                   \
    /* Move the result to the beginning of the dest vec */                                              \
    /* Currently, dest_ind points to one before the start of the legit destination values */            \
    if (dest_ind >= 0 && !status)                                                                       \
    {                                                                                                   \
        new_nnz = (new_nnz-1) - dest_ind;                                                               \
        dest_ind++;                                                                                     \
        for (i = 0; i < new_nnz; i++)                                                                   \
            GR_SPV_SWAP_INDS(DEST_VEC, i, dest_ind + i, CTX);                                           \
    }                                                                                                   \
    (DEST_VEC)->nnz = new_nnz;                                                                          \
    return status;

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_neg_other(gr_ptr res, gr_srcptr src, gr_ctx_t src_ctx, gr_ctx_t ctx) { return (gr_set_other(res, src, src_ctx, ctx) | gr_neg(res, res, ctx)); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_negmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return (gr_mul(res, x, y, ctx) | gr_neg(res, res, ctx)); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_negmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return (gr_mul_si(res, x, y, ctx) | gr_neg(res, res, ctx)); }

#define GR_SPV_RFL_UOP(F, Y, Y_ind) F(GR_ENTRY(res->entries, dest_ptr, sz), GR_ENTRY(Y->entries, Y_ind, sz), ctx)
#define GR_SPV_RFL_BOP(F, Y, Y_ind, Z, Z_ind) F(GR_ENTRY(res->entries, dest_ptr, sz), GR_ENTRY(Y->entries, Y_ind, sz), GR_ENTRY(Z->entries, Z_ind, sz), ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_update(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP(gr_set, src, b_ind), GR_SPV_RFL_UOP(gr_set, src, b_ind), res, res, src, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, src1, a_ind), GR_SPV_RFL_UOP(gr_set, src2, b_ind), GR_SPV_RFL_BOP(gr_add, src1, a_ind, src2, b_ind), res, src1, src2, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, src1, a_ind), GR_SPV_RFL_UOP(gr_neg, src2, b_ind), GR_SPV_RFL_BOP(gr_sub, src1, a_ind, src2, b_ind), res, src1, src2, ctx); }

#define GR_SPV_RFL_UOP_OTHER(F, Y, Y_ind, CTX2) F(GR_ENTRY(res->entries, dest_ptr, sz), GR_ENTRY(Y->entries, Y_ind, CTX2->sizeof_elem), CTX2, ctx)
#define GR_SPV_RFL_BOP_OTHER(F, Y, Y_ind, Z, Z_ind, CTX2) F(GR_ENTRY(res->entries, dest_ptr, sz), GR_ENTRY(Y->entries, Y_ind, sz), GR_ENTRY(Z->entries, Z_ind, CTX2->sizeof_elem), CTX2, ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_add_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, src1, a_ind), GR_SPV_RFL_UOP_OTHER(gr_set, src2, b_ind, ctx2), GR_SPV_RFL_BOP_OTHER(gr_add, src1, a_ind, src2, b_ind, ctx2), res, src1, src2, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sub_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, src1, a_ind), GR_SPV_RFL_UOP_OTHER(gr_neg_other, src2, a_ind, ctx2), GR_SPV_RFL_BOP_OTHER(gr_sub, src1, a_ind, src2, b_ind, ctx2), res, src1, src2, ctx); }

#define GR_SPV_RFL_OTHER_BOP(F, Y, Y_ind, CTX2, Z, Z_ind) F(GR_ENTRY(res->entries, dest_ptr, sz), GR_ENTRY(Y->entries, Y_ind, CTX2->sizeof_elem), GR_ENTRY(Z->entries, Z_ind, sz), ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_add_sparse_vec(gr_sparse_vec_t res, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP_OTHER(gr_set_other, src1, a_ind, ctx1), GR_SPV_RFL_UOP(gr_set, src2, b_ind), GR_SPV_RFL_OTHER_BOP(gr_other_add, src1, a_ind, ctx1, src2, b_ind), res, src1, src2, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_other_sub_sparse_vec(gr_sparse_vec_t res, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP_OTHER(gr_set_other, src1, a_ind, ctx1), GR_SPV_RFL_UOP(gr_neg, src2, b_ind), GR_SPV_RFL_OTHER_BOP(gr_other_sub, src1, a_ind, ctx1, src2, b_ind), res, src1, src2, ctx); }

#define GR_SPV_RFL_UOP_SCALAR(F, Y, Y_ind) F(GR_ENTRY(res->entries, dest_ptr, sz), GR_ENTRY(Y->entries, Y_ind, sz), c, ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP_SCALAR(gr_mul, src, b_ind), GR_SPV_RFL_UOP_SCALAR(gr_addmul, src, b_ind), res, res, src, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP_SCALAR(gr_negmul, src, b_ind), GR_SPV_RFL_UOP_SCALAR(gr_submul, src, b_ind), res, res, src, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_addmul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP_SCALAR(gr_mul_si, src, b_ind), GR_SPV_RFL_UOP_SCALAR(gr_addmul_si, src, b_ind), res, res, src, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_submul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPV_RFL_TEMPLATE(GR_SPV_RFL_UOP(gr_set, res, a_ind), GR_SPV_RFL_UOP_SCALAR(gr_negmul_si, src, b_ind), GR_SPV_RFL_UOP_SCALAR(gr_submul_si, src, b_ind), res, res, src, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_neg(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx) { return (gr_sparse_vec_set(res, src, ctx) | _gr_vec_neg(res->entries, res->entries, res->nnz, ctx)); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_sum(gr_ptr res, const gr_sparse_vec_t vec, gr_ctx_t ctx) { return _gr_vec_sum(res, vec->entries, vec->nnz, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_nz_product(gr_ptr res, const gr_sparse_vec_t vec, gr_ctx_t ctx) { return _gr_vec_product(res, vec->entries, vec->nnz, ctx); }



/***** TODO: These plus dense, dense, sparse functions 
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_MUL)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_DIV)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_DIVEXACT)(res, src1, src2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_VEC_OP(ctx, VEC_POW)(res, src1, src2, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_OP(ctx, VEC_DOT)(res, initial, subtract, vec1, vec2, len, ctx); }

*/




#define GR_SPARSE_VEC_DENSE_VEC_OP(dense_vec_op, res, src, c, ctx) \
{  \
    if(res->length != src->length) \
    { \
        return GR_DOMAIN; \
    } \
    if(res != src) { \
    { \
        gr_sparse_vec_fit_nnz(res, src->nnz, ctx); \
        res->nnz = src->nnz; \
        memcopy(res->cols, src->cols, nnz*sizeof(slong)); \

    } \
    return dense_vec_op(res->entries, src->entries, src->nnz, c, ctx); \
}

int gr_sparse_vec_mul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar, res, src, c, ctx) } 
int gr_sparse_vec_mul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_si, res, src, c, ctx) }
int gr_sparse_vec_mul_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_ui, res, src, c, ctx) }
int gr_sparse_vec_mul_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpz, res, src, c, ctx) }
int gr_sparse_vec_mul_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpq, res, src, c, ctx) }
int gr_sparse_vec_mul_scalar_2exp_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_2exp_si, res, src, c, ctx) }
int gr_sparse_vec_div_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar, res, src, c, ctx) } 
int gr_sparse_vec_div_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_si, res, src, c, ctx) }
int gr_sparse_vec_div_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_ui, res, src, c, ctx) }
int gr_sparse_vec_div_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_fmpz, res, src, c, ctx) }
int gr_sparse_vec_div_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_fmpq, res, src, c, ctx) }
int gr_sparse_vec_divexact_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar, res, src, c, ctx) } 
int gr_sparse_vec_divexact_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_si, res, src, c, ctx) }
int gr_sparse_vec_divexact_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_ui, res, src, c, ctx) }
int gr_sparse_vec_divexact_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpz, res, src, c, ctx) }
int gr_sparse_vec_divexact_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq c, gr_ctx_t ctx) { GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpq, res, src, c, ctx) }



/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/
/******************************************************************/


/*

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

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_MUL_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_div_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_DIV_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_divexact_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_DIVEXACT_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_pow_scalar_other(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) { return GR_SPARSE_VEC_OP_SCALAR_OTHER(ctx, VEC_POW_SCALAR_OTHER)(vec1, vec2, len, c, cctx, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_mul_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_MUL_VEC)(vec1, c, cctx, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_div_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_DIV_VEC)(vec1, c, cctx, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_divexact_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_DIVEXACT_VEC)(vec1, c, cctx, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_scalar_other_pow_vec(gr_sparse_vec_t vec1, gr_srcptr c, gr_ctx_t cctx, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SCALAR_OTHER_OP_VEC(ctx, SCALAR_OTHER_POW_VEC)(vec1, c, cctx, vec2, len, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_mul_scalar_2exp_si(gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, slong c, gr_ctx_t ctx) { return GR_SPARSE_VEC_SCALAR_OP_SI(ctx, VEC_MUL_SCALAR_2EXP_SI)(vec1, vec2, len, c, ctx); }

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

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_OP(ctx, VEC_DOT)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_OP(ctx, VEC_DOT_REV)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const slong * vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_SI_OP(ctx, VEC_DOT_SI)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const ulong * vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_UI_OP(ctx, VEC_DOT_UI)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int gr_sparse_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const fmpz * vec2, slong len, gr_ctx_t ctx) { return GR_SPARSE_VEC_DOT_FMPZ_OP(ctx, VEC_DOT_FMPZ)(res, initial, subtract, vec1, vec2, len, ctx); }
*/



#ifdef __cplusplus
}
#endif

#endif
