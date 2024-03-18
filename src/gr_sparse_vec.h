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
    ulong *inds;
    slong length;
    slong nnz;
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
    _gr_vec_clear(vec->entries, vec->alloc, ctx);
    flint_free(vec->inds);
    flint_free(vec->entries);
    memset(vec, 0, sizeof(gr_sparse_vec_t));
}

#define GR_SPARSE_VEC_IND(vec, nzidx) (vec)->inds[nzidx]
#define GR_SPARSE_VEC_ENTRY(vec, nzidx, sz) GR_ENTRY((vec)->entries, nzidx, sz)

GR_SPARSE_VEC_INLINE ulong * gr_sparse_vec_ind_ptr(gr_sparse_vec_t vec, slong nzidx, gr_ctx_t ctx) { return vec->inds + nzidx; }
GR_SPARSE_VEC_INLINE const ulong * gr_sparse_vec_ind_srcptr(const gr_sparse_vec_t vec, slong nzidx, gr_ctx_t ctx) { return vec->inds + nzidx; }
GR_SPARSE_VEC_INLINE gr_ptr gr_sparse_vec_entry_ptr(gr_sparse_vec_t vec, slong nzidx, gr_ctx_t ctx) { return GR_SPARSE_VEC_ENTRY(vec, nzidx, ctx->sizeof_elem); }
GR_SPARSE_VEC_INLINE gr_srcptr gr_sparse_vec_entry_srcptr(const gr_sparse_vec_t vec, slong nzidx, gr_ctx_t ctx) { return GR_SPARSE_VEC_ENTRY(vec, nzidx, ctx->sizeof_elem); }

void gr_sparse_vec_fit_nnz(gr_sparse_vec_t vec, slong nnz, gr_ctx_t ctx);
void gr_sparse_vec_shrink_to_nnz(gr_sparse_vec_t vec, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE slong gr_sparse_vec_length(const gr_sparse_vec_t vec) { return vec->length; }
void gr_sparse_vec_set_length(gr_sparse_vec_t vec, slong len, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE slong gr_sparse_vec_nnz(const gr_sparse_vec_t vec) { return vec->nnz; }

WARN_UNUSED_RESULT int gr_sparse_vec_set(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_set_from_entries(gr_sparse_vec_t vec, ulong * inds, gr_srcptr entries, slong nnz, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_set_from_entries_sorted_deduped(gr_sparse_vec_t vec, ulong * sorted_deduped_inds, gr_srcptr entries, slong nnz, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE void gr_sparse_vec_swap(gr_sparse_vec_t res, gr_sparse_vec_t src, gr_ctx_t ctx)
{
    FLINT_SWAP(slong, res->alloc, src->alloc);
    FLINT_SWAP(gr_ptr, res->entries, src->entries);
    FLINT_SWAP(ulong *, res->inds, src->inds);
    FLINT_SWAP(slong, res->length, src->length);
    FLINT_SWAP(slong, res->nnz, src->nnz);
    
}
GR_SPARSE_VEC_INLINE void gr_sparse_vec_zero(gr_sparse_vec_t vec, gr_ctx_t ctx) { vec->nnz = 0; }
WARN_UNUSED_RESULT int gr_sparse_vec_set_vec(gr_sparse_vec_t vec, gr_srcptr src, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_vec_set_sparse_vec(gr_ptr vec, gr_sparse_vec_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_slice(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong ind_start, slong ind_end, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_permute_inds(gr_sparse_vec_t vec, const gr_sparse_vec_t src, slong * p, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_find_entry(gr_ptr res, gr_sparse_vec_t vec, slong ind, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_sparse_vec_set_entry(gr_sparse_vec_t vec, slong ind, gr_srcptr entry, gr_ctx_t ctx);

truth_t gr_sparse_vec_equal(const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx);
GR_SPARSE_VEC_INLINE truth_t gr_sparse_vec_is_zero(const gr_sparse_vec_t vec, gr_ctx_t ctx) { return (vec->nnz == 0 ? T_TRUE : T_FALSE); }

int gr_sparse_vec_write_nz(gr_stream_t out, const gr_sparse_vec_t vec, gr_ctx_t ctx);
int gr_sparse_vec_print_nz(const gr_sparse_vec_t vec, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_sparse_vec_randtest(gr_sparse_vec_t vec, double density, slong len, flint_rand_t state, gr_ctx_t ctx);

slong _gr_sparse_vec_count_unique_inds(const ulong *inds0, slong nnz0, const ulong *inds1, slong nnz1);


#define GR_SPV_SWAP_INDS(VEC, I, J, SZ, CTX) \
{                                            \
    slong _temp = (VEC)->inds[I];            \
    (VEC)->inds[I] = (VEC)->inds[J];         \
    (VEC)->inds[J] = _temp;                  \
    gr_swap(GR_ENTRY((VEC)->entries, (I), (SZ)), GR_ENTRY((VEC)->entries, (J), (SZ)), (CTX)); \
}

/*
    Binary arithmetic operations on sparse vectors have a common pattern of
    riffling through the two input vectors A_VEC and B_VEC to produce an
    output vector DEST_VEC with indices in sorted order. As the iteration
    proceeds, for a given index, we need one function FUNC_A if only A_VEC
    has a nonzero element at that index, another function FUNC_B if only
    B_VEC has a nonzero element, and a third function FUNC_AB if both are nonzero.
    During this process, we maintain a triple of indices a_nzidx, b_nzidx, and
    dest_ind, which are used when calling this macro to indicate what element(s)
    to refer to when calling the appropriate function.
*/
#define GR_SPV_RFL_TEMPLATE(FUNC_A, FUNC_B, FUNC_AB, DEST_VEC, A_VEC, B_VEC, CTX)                       \
    int status;                                                                                         \
    slong sz, new_nnz, a_nzidx, b_nzidx, dest_ind, a_nnz, b_nnz, i;                                         \
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
    a_nzidx = a_nnz-1;                                                                                    \
    b_nzidx = b_nnz-1;                                                                                    \
    dest_ind = new_nnz-1;                                                                               \
    while (a_nzidx >= 0 && b_nzidx >= 0 && status == GR_SUCCESS)                                            \
    {                                                                                                   \
        a_ind = (A_VEC)->inds[a_nzidx];                                                                   \
        b_ind = (B_VEC)->inds[b_nzidx];                                                                   \
        if (a_ind > b_ind)                                                                              \
        {                                                                                               \
            status |= (FUNC_A);                                                                         \
            (DEST_VEC)->inds[dest_ind] = a_ind;                                                         \
            a_nzidx--;                                                                                    \
        }                                                                                               \
        else if (b_ind > a_ind)                                                                         \
        {                                                                                               \
            status |= (FUNC_B);                                                                         \
            (DEST_VEC)->inds[dest_ind] = b_ind;                                                         \
            b_nzidx--;                                                                                    \
        }                                                                                               \
        else                                                                                            \
        {                                                                                               \
            status |= (FUNC_AB);                                                                        \
            (DEST_VEC)->inds[dest_ind] = a_ind;                                                         \
            a_nzidx--;                                                                                    \
            b_nzidx--;                                                                                    \
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

/* We need some convenience functions for certain simple operations. */
GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_neg_other(gr_ptr res, gr_srcptr src, gr_ctx_t src_ctx, gr_ctx_t ctx)
{ return (gr_set_other(res, src, src_ctx, ctx) | gr_neg(res, res, ctx)); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_negmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{ return (gr_mul(res, x, y, ctx) | gr_neg(res, res, ctx)); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_negmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{ return (gr_mul_si(res, x, y, ctx) | gr_neg(res, res, ctx)); }

// Convenience macros for applying a unary or binary function
#define GR_SPV_RFL_UOP(F, Y, Y_ind) \
    F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, sz), ctx)
#define GR_SPV_RFL_BOP(F, Y, Y_ind, Z, Z_ind) \
    F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, sz), GR_ENTRY(Z->entries, Z_ind, sz), ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_update(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx)
{
    // All function are set to src value if possible, otherwise res value
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, res, a_nzidx), 
        GR_SPV_RFL_UOP(gr_set, src, b_nzidx), 
        GR_SPV_RFL_UOP(gr_set, src, b_nzidx),
        res, res, src, ctx
    ); 
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_add(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, src1, a_nzidx), 
        GR_SPV_RFL_UOP(gr_set, src2, b_nzidx), 
        GR_SPV_RFL_BOP(gr_add, src1, a_nzidx, src2, b_nzidx), 
        res, src1, src2, ctx
    );
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_sub(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx) 
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, src1, a_nzidx),
        GR_SPV_RFL_UOP(gr_neg, src2, b_nzidx),
        GR_SPV_RFL_BOP(gr_sub, src1, a_nzidx, src2, b_nzidx),
        res, src1, src2, ctx
    );
}

// Need a unary function to assign zero
#define GR_SPV_RFL_ZERO gr_zero(GR_ENTRY(res->entries, dest_ind, sz), ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_mul(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, slong len, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_BOP(gr_mul, src1, a_nzidx, src2, b_nzidx),
        res, src1, src2, ctx
    );
}

// Analogous macros for applying a unary or binary function between two contexts
#define GR_SPV_RFL_UOP_OTHER(F, Y, Y_ind, CTX2) \
    F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, CTX2->sizeof_elem), CTX2, ctx)
#define GR_SPV_RFL_BOP_OTHER(F, Y, Y_ind, Z, Z_ind, CTX2) \
    F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, sz), GR_ENTRY(Z->entries, Z_ind, CTX2->sizeof_elem), CTX2, ctx)
#define GR_SPV_RFL_OTHER_BOP(F, Y, Y_ind, CTX2, Z, Z_ind) \
    F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, CTX2->sizeof_elem), (CTX2), GR_ENTRY(Z->entries, Z_ind, sz), ctx)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_add_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, src1, a_nzidx), 
        GR_SPV_RFL_UOP_OTHER(gr_set_other, src2, b_nzidx, ctx2),
        GR_SPV_RFL_BOP_OTHER(gr_add_other, src1, a_nzidx, src2, b_nzidx, ctx2),
        res, src1, src2, ctx
    );
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_sub_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP(gr_set, src1, a_nzidx),
        GR_SPV_RFL_UOP_OTHER(gr_neg_other, src2, a_nzidx, ctx2),
        GR_SPV_RFL_BOP_OTHER(gr_sub_other, src1, a_nzidx, src2, b_nzidx, ctx2),
        res, src1, src2, ctx
    );
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_mul_other(gr_sparse_vec_t res, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_BOP_OTHER(gr_mul_other, src1, a_nzidx, src2, b_nzidx, ctx2),
        res, src1, src2, ctx
    );
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_other_add_sparse_vec(gr_sparse_vec_t res, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP_OTHER(gr_set_other, src1, a_nzidx, ctx1),
        GR_SPV_RFL_UOP(gr_set, src2, b_nzidx),
        GR_SPV_RFL_OTHER_BOP(gr_other_add, src1, a_nzidx, ctx1, src2, b_nzidx),
        res, src1, src2, ctx
    );
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_other_sub_sparse_vec(gr_sparse_vec_t res, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_UOP_OTHER(gr_set_other, src1, a_nzidx, ctx1),
        GR_SPV_RFL_UOP(gr_neg, src2, b_nzidx),
        GR_SPV_RFL_OTHER_BOP(gr_other_sub, src1, a_nzidx, ctx1, src2, b_nzidx),
        res, src1, src2, ctx
    );
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_other_mul_sparse_vec(gr_sparse_vec_t res, const gr_sparse_vec_t src1, gr_ctx_t ctx1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    GR_SPV_RFL_TEMPLATE(
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_ZERO,
        GR_SPV_RFL_OTHER_BOP(gr_other_mul, src1, a_nzidx, ctx1, src2, b_nzidx),
        res, src1, src2, ctx
    );
}

#define GR_SPV_RFL_UOP_SCALAR(F, Y, Y_ind, C, CTX) \
    F(GR_ENTRY(res->entries, dest_ind, sz), GR_ENTRY(Y->entries, Y_ind, sz), C, CTX)

/*
    Subtemplate for doing accumulated operation from one sparse vector into another.
    Used to do res += c * src and res -= c * src, for c a scalar.
*/
#define GR_SPV_ACCUM_TEMPLATE(ELEM_OP, ELEM_ACCUM_OP, RES, SRC, C, CTX) \
    GR_SPV_RFL_TEMPLATE(                                                \
        GR_SPV_RFL_UOP(gr_set, RES, a_nzidx),                             \
        GR_SPV_RFL_UOP_SCALAR(ELEM_OP, SRC, b_nzidx, C, CTX),             \
        GR_SPV_RFL_UOP_SCALAR(ELEM_ACCUM_OP, SRC, b_nzidx, C, CTX),       \
        RES, RES, SRC, CTX\
    )

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_addmul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_TEMPLATE(gr_mul, gr_addmul, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_submul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_TEMPLATE(gr_negmul, gr_submul, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_addmul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_TEMPLATE(gr_mul_si, gr_addmul_si, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_submul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_TEMPLATE(gr_negmul_si, gr_submul_si, res, src, c, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_neg(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_ctx_t ctx)
{ return (gr_sparse_vec_set(res, src, ctx) | _gr_vec_neg(res->entries, res->entries, res->nnz, ctx)); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_sum(gr_ptr res, const gr_sparse_vec_t vec, gr_ctx_t ctx)
{ return _gr_vec_sum(res, vec->entries, vec->nnz, ctx); }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_sparse_vec_nz_product(gr_ptr res, const gr_sparse_vec_t vec, gr_ctx_t ctx)
{ return _gr_vec_product(res, vec->entries, vec->nnz, ctx); }


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

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_update_sparse_vec_nz(gr_ptr dres, const gr_sparse_vec_t src, gr_ctx_t ctx)
{ GR_SPV_INTO_DENSE_TEMPLATE(gr_set(GR_ENTRY(dres, src->inds[i], sz), GR_ENTRY(src->entries, i, sz), ctx), src, ctx) }

/* Sub-macro for applying operation (dense, sparse) -> dense */

#define GR_SPV_OP_ON_DENSE_TEMPLATE(ELEM_OP, DRES, DVEC, SVEC, CTX) \
    GR_SPV_INTO_DENSE_TEMPLATE(ELEM_OP(GR_ENTRY((DRES), (SVEC)->inds[i], sz), GR_ENTRY(DVEC, (SVEC)->inds[i], sz), GR_ENTRY((SVEC)->entries, i, sz), CTX), SVEC, CTX)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_add_sparse_vec(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
{ GR_SPV_OP_ON_DENSE_TEMPLATE(gr_add, dres, dvec1, svec2, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_sub_sparse_vec(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
{ GR_SPV_OP_ON_DENSE_TEMPLATE(gr_sub, dres, dvec1, svec2, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_mul_sparse_vec_nz(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
{ GR_SPV_OP_ON_DENSE_TEMPLATE(gr_mul, dres, dvec1, svec2, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_div_sparse_vec_nz(gr_ptr dres, gr_srcptr dvec1, const gr_sparse_vec_t svec2, gr_ctx_t ctx)
{ GR_SPV_OP_ON_DENSE_TEMPLATE(gr_div, dres, dvec1, svec2, ctx) }

/* Sub-macro for accumulating operation on dense from sparse */

#define GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(ELEM_OP, DRES, SVEC, C, CTX) \
    GR_SPV_INTO_DENSE_TEMPLATE(ELEM_OP(GR_ENTRY(DRES, (SVEC)->inds[i], sz), GR_ENTRY((SVEC)->entries, i, sz), C, CTX), SVEC, CTX)

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_addmul_sparse_vec_scalar(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_addmul, dres, svec, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_submul_sparse_vec_scalar(gr_ptr dres, const gr_sparse_vec_t svec, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_submul, dres, svec, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_addmul_sparse_vec_scalar_si(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_addmul_si, dres, svec, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_submul_sparse_vec_scalar_si(gr_ptr dres, const gr_sparse_vec_t svec, slong c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_submul_si, dres, svec, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_addmul_sparse_vec_scalar_fmpz(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_addmul_fmpz, dres, svec, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int 
gr_vec_submul_sparse_vec_scalar_fmpz(gr_ptr dres, const gr_sparse_vec_t svec, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPV_ACCUM_INTO_DENSE_TEMPLATE(gr_submul_fmpz, dres, svec, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const gr_sparse_vec_t vec2, gr_ctx_t ctx)
{
    int status;
    slong nzidx1, nzidx2;
    slong sz = ctx->sizeof_elem;
    
    if (vec1->length != vec2->length)
    {
        return GR_DOMAIN;
    }
    status = gr_set(res, initial, ctx);
    for (nzidx1 = 0, nzidx2 = 0; nzidx1 < vec1->nnz && nzidx2 < vec2->nnz; )
    {
        if (vec1->inds[nzidx1] < vec2->inds[nzidx2])
        {
            nzidx1++;
        }
        else if (vec1->inds[nzidx1] > vec2->inds[nzidx2])
        {
            nzidx2++;
        }
        else if (subtract)
        {
            status |= gr_submul(res, GR_ENTRY(vec1->entries, nzidx1, sz), GR_ENTRY(vec2->entries, nzidx1, sz), ctx);
        }
        else
        {
            status |= gr_addmul(res, GR_ENTRY(vec1->entries, nzidx1, sz), GR_ENTRY(vec2->entries, nzidx1, sz), ctx);
        }
    }
    return status;
}

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_dot_vec(gr_ptr res, gr_srcptr initial, int subtract, const gr_sparse_vec_t vec1, const gr_vec_t vec2, gr_ctx_t ctx)
{
    int status;
    slong nzidx;
    slong sz = ctx->sizeof_elem;
    
    if (vec1->length != vec2->length)
    {
        return GR_DOMAIN;
    }
    status = gr_set(res, initial, ctx);
    for (nzidx = 0; nzidx < vec1->nnz; nzidx++)
    {
        if (subtract)
        {
            status |= gr_submul(res, GR_ENTRY(vec1->entries, nzidx, sz), GR_ENTRY(vec2->entries, vec1->inds[nzidx], sz), ctx);
        }
        else
        {
            status |= gr_addmul(res, GR_ENTRY(vec1->entries, nzidx, sz), GR_ENTRY(vec2->entries, vec1->inds[nzidx], sz), ctx);
        }
    }
    return status;
}

#define GR_SPARSE_VEC_DENSE_VEC_OP(dense_vec_op, res, src, c, ctx)     \
    if(res->length != src->length)                                     \
    {                                                                  \
        return GR_DOMAIN;                                              \
    }                                                                  \
    if(res != src)                                                     \
    {                                                                  \
        gr_sparse_vec_fit_nnz(res, src->nnz, ctx);                     \
        res->nnz = src->nnz;                                           \
        memcpy(res->inds, src->inds, src->nnz*sizeof(slong));          \
    }                                                                  \
    return dense_vec_op(res->entries, src->entries, src->nnz, c, ctx); \

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_mul_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar, res, src, c, ctx) } 

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_mul_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_si, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_mul_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_ui, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_mul_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpz, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_mul_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpq, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_mul_scalar_2exp_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_mul_scalar_2exp_si, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_div_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar, res, src, c, ctx) } 

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_div_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_si, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_div_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_ui, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_div_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_fmpz, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_div_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_div_scalar_fmpq, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_divexact_scalar(gr_sparse_vec_t res, const gr_sparse_vec_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar, res, src, c, ctx) } 

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_divexact_scalar_si(gr_sparse_vec_t res, const gr_sparse_vec_t src, slong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_si, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_divexact_scalar_ui(gr_sparse_vec_t res, const gr_sparse_vec_t src, ulong c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_ui, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_divexact_scalar_fmpz(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpz, res, src, c, ctx) }

GR_SPARSE_VEC_INLINE WARN_UNUSED_RESULT int
gr_sparse_vec_divexact_scalar_fmpq(gr_sparse_vec_t res, const gr_sparse_vec_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_SPARSE_VEC_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpq, res, src, c, ctx) }



#ifdef __cplusplus
}
#endif

#endif
