/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "gr_sparse_vec.h"

/*
    Binary arithmetic operations on sparse vectors have a common pattern of
    riffling through the two input vectors A_VEC and B_VEC to produce an
    output vector DEST_VEC with indices in sorted order. As the iteration
    proceeds, for a given index, we need one function FUNC_A if only A_VEC
    has a nonzero element at that index, another function FUNC_B if only
    B_VEC has a nonzero element, and a third function FUNC_AB if both are nonzero.
    During this process, we maintain a triple of indices a_nz_idx, b_nz_idx, and
    dest_nz_idx, which are used when calling this macro to indicate what element(s)
    to refer to when calling the appropriate function.
*/

// Sub macro to swap two entries with their associated indices
#define GR_SPV_SWAP_INDS(VEC, I, J, SZ, CTX) \
{                                            \
    slong _temp = (VEC)->inds[I];            \
    (VEC)->inds[I] = (VEC)->inds[J];         \
    (VEC)->inds[J] = _temp;                  \
    gr_swap(GR_ENTRY((VEC)->nzs, (I), (SZ)), GR_ENTRY((VEC)->nzs, (J), (SZ)), (CTX)); \
}

#define GR_SPV_RFL_TEMPLATE(FUNC_A, FUNC_B, FUNC_AB, DEST_VEC, A_VEC, B_VEC, CTX)                       \
    int status;                                                                                         \
    slong sz, new_nnz, a_nz_idx, b_nz_idx, dest_nz_idx, a_nnz, b_nnz, i;                                         \
    ulong a_ind, b_ind;                                                                                 \
    if ((DEST_VEC)->length != (A_VEC)->length || (A_VEC)->length != (B_VEC)->length)                    \
        return GR_DOMAIN;                                                                               \
    status = GR_SUCCESS;                                                                                \
    sz = (CTX)->sizeof_elem;                                                                            \
    a_nnz = (A_VEC)->nnz;                                                                               \
    b_nnz = (B_VEC)->nnz;                                                                               \
    new_nnz = _gr_sparse_vec_count_unique_inds((A_VEC)->inds, a_nnz, (B_VEC)->inds, b_nnz);             \
    gr_sparse_vec_fit_nnz((DEST_VEC), new_nnz, (CTX));                                                  \
    /* We go backward through the destination, because it might be an in-place operation on a source */ \
    a_nz_idx = a_nnz-1;                                                                                    \
    b_nz_idx = b_nnz-1;                                                                                    \
    dest_nz_idx = new_nnz-1;                                                                               \
    while (a_nz_idx >= -1 && b_nz_idx >= -1 && status == GR_SUCCESS)                                       \
    {                                                                                                      \
        if (a_nz_idx == -1 && b_nz_idx == -1) break;                                                       \
        a_ind = (A_VEC)->inds[a_nz_idx];                                                                   \
        b_ind = (B_VEC)->inds[b_nz_idx];                                                                   \
        if (b_nz_idx == -1 || (a_nz_idx >= 0 && a_ind > b_ind))                                                                              \
        {                                                                                               \
            status |= (FUNC_A);                                                                         \
            (DEST_VEC)->inds[dest_nz_idx] = a_ind;                                                         \
            a_nz_idx--;                                                                                    \
        }                                                                                               \
        else if (a_nz_idx == -1 || b_ind > a_ind)                                                                         \
        {                                                                                               \
            status |= (FUNC_B);                                                                         \
            (DEST_VEC)->inds[dest_nz_idx] = b_ind;                                                         \
            b_nz_idx--;                                                                                    \
        }                                                                                               \
        else                                                                                            \
        {                                                                                               \
            status |= (FUNC_AB);                                                                        \
            (DEST_VEC)->inds[dest_nz_idx] = a_ind;                                                         \
            a_nz_idx--;                                                                                    \
            b_nz_idx--;                                                                                    \
        }                                                                                               \
        if (T_TRUE != gr_is_zero(GR_ENTRY((DEST_VEC)->nzs, dest_nz_idx, sz), (CTX)))                   \
            dest_nz_idx--;                                                                                 \
    }                                                                                                   \
    /* Move the result to the beginning of the dest vec */                                              \
    /* Currently, dest_nz_idx points to one before the start of the legit destination values */            \
    if (dest_nz_idx >= 0 && !status)                                                                       \
    {                                                                                                   \
        new_nnz = (new_nnz-1) - dest_nz_idx;                                                               \
        dest_nz_idx++;                                                                                     \
        for (i = 0; i < new_nnz; i++)                                                                   \
            GR_SPV_SWAP_INDS(DEST_VEC, i, dest_nz_idx + i, sz, CTX);                                       \
    }                                                                                                   \
    (DEST_VEC)->nnz = new_nnz;                                                                          \
    return status;

/* We need some convenience functions for certain simple operations. */
int gr_neg_other(gr_ptr dst, gr_srcptr src, gr_ctx_t src_ctx, gr_ctx_t ctx)
{ return (gr_set_other(dst, src, src_ctx, ctx) | gr_neg(dst, dst, ctx)); }

int gr_negmul(gr_ptr dst, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{ return (gr_mul(dst, x, y, ctx) | gr_neg(dst, dst, ctx)); }

int gr_negmul_si(gr_ptr dst, gr_srcptr x, slong y, gr_ctx_t ctx)
{ return (gr_mul_si(dst, x, y, ctx) | gr_neg(dst, dst, ctx)); }

// Need a zero-ary function to assign zero
#define GR_SPV_RFL_ZERO gr_zero(GR_ENTRY(dst->nzs, dest_nz_idx, sz), ctx)

// Convenience macros for applying a unary or binary function
#define GR_SPV_RFL_UOP(F, DST, DST_IND, SRC, SRC_IND) \
    F(GR_ENTRY((DST)->nzs, DST_IND, sz), GR_ENTRY((SRC)->nzs, SRC_IND, sz), ctx)
#define GR_SPV_RFL_BOP(F, DST, DST_IND, SRC1, SRC1_IND, SRC2, SRC2_IND) \
    F(GR_ENTRY((DST)->nzs, DST_IND, sz), GR_ENTRY((SRC1)->nzs, SRC1_IND, sz), GR_ENTRY((SRC2)->nzs, SRC2_IND, sz), ctx)
#define GR_SPV_RFL_BOP_SCALAR(F, DST, DST_IND, SRC, SRC_IND, C, CTX) \
    F(GR_ENTRY((DST)->nzs, DST_IND, sz), GR_ENTRY((SRC)->nzs, SRC_IND, sz), C, CTX)

// Analogous macros for applying a unary or binary function between two contexts
#define GR_SPV_RFL_UOP_OTHER(F, DST, DST_IND, SRC, SRC_IND, CTX2) \
    F(GR_ENTRY((DST)->nzs, DST_IND, sz), GR_ENTRY((SRC)->nzs, SRC_IND, (CTX2)->sizeof_elem), CTX2, ctx)
#define GR_SPV_RFL_BOP_OTHER(F, DST, DST_IND, SRC1, SRC1_IND, SRC2, SRC2_IND, CTX2) \
    F(GR_ENTRY((DST)->nzs, DST_IND, sz), GR_ENTRY((SRC1)->nzs, SRC1_IND, sz), GR_ENTRY((SRC2)->nzs, SRC2_IND, (CTX2)->sizeof_elem), CTX2, ctx)
#define GR_SPV_RFL_OTHER_BOP(F, DST, DST_IND, SRC1, SRC1_IND, CTX2, SRC2, SRC2_IND) \
    F(GR_ENTRY((DST)->nzs, DST_IND, sz), GR_ENTRY((SRC1)->nzs, SRC1_IND, CTX2->sizeof_elem), (CTX2), GR_ENTRY((SRC2)->nzs, SRC2_IND, sz), ctx)

/*
    Subtemplate for doing accumulated operation from one sparse vector into another.
    Used to do dst += c * src and dst -= c * src, for c a scalar.
*/
#define GR_SPV_ACCUM_TEMPLATE(ELEM_OP, ELEM_ACCUM_OP, DST, SRC, C, CTX) \
    GR_SPV_RFL_TEMPLATE(                                                \
        GR_SPV_RFL_UOP(gr_set, DST, dest_nz_idx, DST, a_nz_idx),                             \
        GR_SPV_RFL_BOP_SCALAR(ELEM_OP, DST, dest_nz_idx, SRC, b_nz_idx, C, CTX),             \
        GR_SPV_RFL_UOP(gr_set, DST, dest_nz_idx, DST, a_nz_idx) | GR_SPV_RFL_BOP_SCALAR(ELEM_ACCUM_OP, DST, dest_nz_idx, SRC, b_nz_idx, C, CTX),       \
        DST, DST, SRC, CTX\
    )

int gr_sparse_vec_neg(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx)
{ return (gr_sparse_vec_set(dst, src, ctx) | _gr_vec_neg(dst->nzs, dst->nzs, dst->nnz, ctx)); }

int
gr_sparse_vec_update(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, dst, a_nz_idx), 
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, src, b_nz_idx), 
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, src, b_nz_idx),
        dst, dst, src, ctx
    ); 
}

int
gr_sparse_vec_add(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, src1, a_nz_idx), 
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, src2, b_nz_idx), 
        GR_SPV_RFL_BOP(gr_add, dst, dest_nz_idx, src1, a_nz_idx, src2, b_nz_idx), 
        dst, src1, src2, ctx
    );
}

int
gr_sparse_vec_sub(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx) 
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, src1, a_nz_idx),
        GR_SPV_RFL_UOP(gr_neg, dst, dest_nz_idx, src2, b_nz_idx),
        GR_SPV_RFL_BOP(gr_sub, dst, dest_nz_idx, src1, a_nz_idx, src2, b_nz_idx),
        dst, src1, src2, ctx
    );
}

int gr_sparse_vec_mul(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_BOP(gr_mul, dst, dest_nz_idx, src1, a_nz_idx, src2, b_nz_idx),
        dst, src1, src2, ctx
    );
}

int
gr_sparse_vec_add_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, src1, a_nz_idx), 
        GR_SPV_RFL_UOP_OTHER(gr_set_other, dst, dest_nz_idx, src2, b_nz_idx, ctx2),
        GR_SPV_RFL_BOP_OTHER(gr_add_other, dst, dest_nz_idx, src1, a_nz_idx, src2, b_nz_idx, ctx2),
        dst, src1, src2, ctx
    );
}

int
gr_sparse_vec_sub_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, src1, a_nz_idx),
        GR_SPV_RFL_UOP_OTHER(gr_neg_other, dst, dest_nz_idx, src2, b_nz_idx, ctx2),
        GR_SPV_RFL_BOP_OTHER(gr_sub_other, dst, dest_nz_idx, src1, a_nz_idx, src2, b_nz_idx, ctx2),
        dst, src1, src2, ctx
    );
}

int
gr_sparse_vec_mul_other(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_BOP_OTHER(gr_mul_other, dst, dest_nz_idx, src1, a_nz_idx, src2, b_nz_idx, ctx2),
        dst, src1, src2, ctx
    );
}

int gr_other_add_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP_OTHER(gr_set_other, dst, dest_nz_idx, src1, a_nz_idx, ctx1),
        GR_SPV_RFL_UOP(gr_set, dst, dest_nz_idx, src2, b_nz_idx),
        GR_SPV_RFL_OTHER_BOP(gr_other_add, dst, dest_nz_idx, src1, a_nz_idx, ctx1, src2, b_nz_idx),
        dst, src1, src2, ctx
    );
}

int gr_other_sub_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP_OTHER(gr_set_other, dst, dest_nz_idx, src1, a_nz_idx, ctx1),
        GR_SPV_RFL_UOP(gr_neg, dst, dest_nz_idx, src2, b_nz_idx),
        GR_SPV_RFL_OTHER_BOP(gr_other_sub, dst, dest_nz_idx, src1, a_nz_idx, ctx1, src2, b_nz_idx),
        dst, src1, src2, ctx
    );
}

int gr_other_mul_sparse_vec(gr_sparse_vec_t dst, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_OTHER_BOP(gr_other_mul, dst, dest_nz_idx, src1, a_nz_idx, ctx1, src2, b_nz_idx),
        dst, src1, src2, ctx
    );
}

int gr_sparse_vec_addmul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_TEMPLATE(gr_mul, gr_addmul, dst, src, c, ctx) }

int gr_sparse_vec_submul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_TEMPLATE(gr_negmul, gr_submul, dst, src, c, ctx) }

int gr_sparse_vec_addmul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_TEMPLATE(gr_mul_si, gr_addmul_si, dst, src, c, ctx) }

int gr_sparse_vec_submul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_TEMPLATE(gr_negmul_si, gr_submul_si, dst, src, c, ctx); }

/**
 * Arithmetic into dense vectors
**/

// Internal macro to update a dense vector by iterating over a sparse one
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

// Sub-macro for applying operation (dense, sparse) -> dense
#define GR_SPV_OP_ON_DENSE_TEMPLATE(ELEM_OP, DRES, DVEC, SVEC, CTX) \
    GR_SPV_INTO_DENSE_TEMPLATE(ELEM_OP(\
        GR_ENTRY((DRES), (SVEC)->inds[i], sz), \
        GR_ENTRY(DVEC, (SVEC)->inds[i], sz), \
        GR_ENTRY((SVEC)->nzs, i, sz), \
    CTX), SVEC, CTX)

// Sub-macro for accumulating operation on dense from sparse
#define GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(ELEM_OP, DRES, SVEC, C, CTX) \
    GR_SPV_INTO_DENSE_TEMPLATE(ELEM_OP(\
        GR_ENTRY(DRES, (SVEC)->inds[i], sz), \
        GR_ENTRY((SVEC)->nzs, i, sz),\
    C, CTX), SVEC, CTX)

int gr_vec_update_sparse_vec_nz(gr_ptr dres, const gr_sparse_vec_t src, gr_ctx_t ctx)
{ GR_SPV_INTO_DENSE_TEMPLATE(gr_set(GR_ENTRY(dres, src->inds[i], sz), GR_ENTRY(src->nzs, i, sz), ctx), src, ctx) }

int gr_vec_add_sparse_vec(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
{ GR_SPV_OP_ON_DENSE_TEMPLATE(gr_add, dres, dvec1, svec2, ctx) }

int gr_vec_sub_sparse_vec(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
{ GR_SPV_OP_ON_DENSE_TEMPLATE(gr_sub, dres, dvec1, svec2, ctx) }

int gr_vec_mul_sparse_vec_nz(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
{ GR_SPV_OP_ON_DENSE_TEMPLATE(gr_mul, dres, dvec1, svec2, ctx) }

int gr_vec_div_sparse_vec_nz(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
{ GR_SPV_OP_ON_DENSE_TEMPLATE(gr_div, dres, dvec1, svec2, ctx) }

int gr_vec_addmul_sparse_vec_scalar(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_addmul, dres, svec, c, ctx) }

int gr_vec_submul_sparse_vec_scalar(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_submul, dres, svec, c, ctx) }

int gr_vec_addmul_sparse_vec_scalar_si(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_addmul_si, dres, svec, c, ctx) }

int gr_vec_submul_sparse_vec_scalar_si(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_submul_si, dres, svec, c, ctx) }

int gr_vec_addmul_sparse_vec_scalar_fmpz(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_addmul_fmpz, dres, svec, c, ctx) }

int gr_vec_submul_sparse_vec_scalar_fmpz(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_submul_fmpz, dres, svec, c, ctx) }

/**
 * Scalar multiplication and division
**/

#define GR_SPARSE_VEC_DENSE_VEC_OP(dense_vec_op, dst, src, c, ctx)     \
    if(dst->length != src->length)                                     \
    {                                                                  \
        return GR_DOMAIN;                                              \
    }                                                                  \
    if(dst != src)                                                     \
    {                                                                  \
        gr_sparse_vec_fit_nnz(dst, src->nnz, ctx);                     \
        dst->nnz = src->nnz;                                           \
        memcpy(dst->inds, src->inds, src->nnz*sizeof(slong));          \
    }                                                                  \
    return dense_vec_op(dst->nzs, src->nzs, src->nnz, c, ctx); \

int
gr_sparse_vec_mul_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar, dst, src, c, ctx) } 

int
gr_sparse_vec_mul_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_si, dst, src, c, ctx) }

int
gr_sparse_vec_mul_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_ui, dst, src, c, ctx) }

int
gr_sparse_vec_mul_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpz, dst, src, c, ctx) }

int
gr_sparse_vec_mul_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpq, dst, src, c, ctx) }

int
gr_sparse_vec_mul_scalar_2exp_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_2exp_si, dst, src, c, ctx) }

int
gr_sparse_vec_div_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar, dst, src, c, ctx) } 

int
gr_sparse_vec_div_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_si, dst, src, c, ctx) }

int
gr_sparse_vec_div_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_ui, dst, src, c, ctx) }

int
gr_sparse_vec_div_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_fmpz, dst, src, c, ctx) }

int
gr_sparse_vec_div_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_fmpq, dst, src, c, ctx) }

int
gr_sparse_vec_divexact_scalar(gr_sparse_vec_t dst, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar, dst, src, c, ctx) } 

int
gr_sparse_vec_divexact_scalar_si(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_si, dst, src, c, ctx) }

int
gr_sparse_vec_divexact_scalar_ui(gr_sparse_vec_t dst, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_ui, dst, src, c, ctx) }

int
gr_sparse_vec_divexact_scalar_fmpz(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpz, dst, src, c, ctx) }

int
gr_sparse_vec_divexact_scalar_fmpq(gr_sparse_vec_t dst, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpq, dst, src, c, ctx) }
