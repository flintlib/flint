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
#include "gr_sparse_mat.h"

#define GR_CSR_MAT_DENSE_VEC_OP(dense_vec_op, dst, src, c, ctx)     \
    if(dst->r != src->r || dst->c != src->c)                           \
    {                                                                  \
        return GR_DOMAIN;                                              \
    }                                                                  \
    if(dst != src)                                                     \
    {                                                                  \
        gr_csr_mat_fit_nnz(dst, src->nnz, ctx);                     \
        dst->nnz = src->nnz;                                           \
        memcpy(dst->rows, src->rows, src->r*sizeof(ulong));          \
        memcpy(dst->cols, src->cols, src->nnz*sizeof(ulong));          \
    }                                                                  \
    return dense_vec_op(dst->nzs, src->nzs, src->nnz, c, ctx); \

int gr_csr_mat_mul_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar, dst, src, c, ctx) } 
int gr_csr_mat_mul_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_si, dst, src, c, ctx) }
int gr_csr_mat_mul_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_ui, dst, src, c, ctx) }
int gr_csr_mat_mul_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpz, dst, src, c, ctx) }
int gr_csr_mat_mul_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpq, dst, src, c, ctx) }
int gr_csr_mat_mul_scalar_2exp_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_2exp_si, dst, src, c, ctx) }
int gr_csr_mat_div_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_div_scalar, dst, src, c, ctx) } 
int gr_csr_mat_div_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_div_scalar_si, dst, src, c, ctx) }
int gr_csr_mat_div_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_div_scalar_ui, dst, src, c, ctx) }
int gr_csr_mat_div_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_div_scalar_fmpz, dst, src, c, ctx) }
int gr_csr_mat_div_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_div_scalar_fmpq, dst, src, c, ctx) }
int gr_csr_mat_divexact_scalar(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar, dst, src, c, ctx) } 
int gr_csr_mat_divexact_scalar_si(gr_csr_mat_t dst, const gr_csr_mat_t src, slong c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar_si, dst, src, c, ctx) }
int gr_csr_mat_divexact_scalar_ui(gr_csr_mat_t dst, const gr_csr_mat_t src, ulong c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar_ui, dst, src, c, ctx) }
int gr_csr_mat_divexact_scalar_fmpz(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpz, dst, src, c, ctx) }
int gr_csr_mat_divexact_scalar_fmpq(gr_csr_mat_t dst, const gr_csr_mat_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_CSR_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpq, dst, src, c, ctx) }

#define GR_LIL_MAT_DENSE_VEC_OP(dense_vec_op, dst, src, c, ctx)   {  \
    int status = GR_SUCCESS;                                          \
    int row;                                                           \
    if(dst->r != src->r || dst->c != src->c)                           \
    {                                                                  \
        return GR_DOMAIN;                                              \
    }                                                                  \
    dst->nnz = src->nnz;                                               \
    for (row = 0; row < src->r; ++row)                                  \
    {                                                                  \
        if(dst != src)                                                     \
        {                                                                  \
            gr_sparse_vec_fit_nnz(&dst->rows[row], src->rows[row].nnz, ctx);                     \
            dst->rows[row].nnz = src->rows[row].nnz;                                           \
            memcpy(dst->rows[row].inds, src->rows[row].inds, src->rows[row].nnz*sizeof(ulong));          \
        }                                                                  \
        status |= dense_vec_op(dst->rows[row].nzs, src->rows[row].nzs, src->rows[row].nnz, c, ctx); \
    }                                                                  \
    return status; \
}

int gr_lil_mat_mul_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar, dst, src, c, ctx) } 
int gr_lil_mat_mul_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_si, dst, src, c, ctx) }
int gr_lil_mat_mul_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_ui, dst, src, c, ctx) }
int gr_lil_mat_mul_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpz, dst, src, c, ctx) }
int gr_lil_mat_mul_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_fmpq, dst, src, c, ctx) }
int gr_lil_mat_mul_scalar_2exp_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_mul_scalar_2exp_si, dst, src, c, ctx) }
int gr_lil_mat_div_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_div_scalar, dst, src, c, ctx) } 
int gr_lil_mat_div_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_div_scalar_si, dst, src, c, ctx) }
int gr_lil_mat_div_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_div_scalar_ui, dst, src, c, ctx) }
int gr_lil_mat_div_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_div_scalar_fmpz, dst, src, c, ctx) }
int gr_lil_mat_div_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_div_scalar_fmpq, dst, src, c, ctx) }
int gr_lil_mat_divexact_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar, dst, src, c, ctx) } 
int gr_lil_mat_divexact_scalar_si(gr_lil_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar_si, dst, src, c, ctx) }
int gr_lil_mat_divexact_scalar_ui(gr_lil_mat_t dst, const gr_lil_mat_t src, ulong c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar_ui, dst, src, c, ctx) }
int gr_lil_mat_divexact_scalar_fmpz(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpz, dst, src, c, ctx) }
int gr_lil_mat_divexact_scalar_fmpq(gr_lil_mat_t dst, const gr_lil_mat_t src, const fmpq_t c, gr_ctx_t ctx)
{ GR_LIL_MAT_DENSE_VEC_OP(_gr_vec_divexact_scalar_fmpq, dst, src, c, ctx) }
