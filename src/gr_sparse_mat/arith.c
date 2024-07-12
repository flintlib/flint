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

#define GR_LIL_MAT_BOP(FUNC, RES, MAT1, MAT2, CTX) \
{ \
    int row; \
    int status = GR_SUCCESS; \
                \
    if (gr_mat_is_compatible(RES, MAT1, CTX) == T_FALSE || gr_mat_is_compatible(RES, MAT2, CTX) == T_FALSE) \
    { \
        return GR_DOMAIN; \
    } \
    (RES)->nnz = 0; \
    for (row = 0; row < (RES)->r; row++) \
{ \
        status |= FUNC(&(RES)->rows[row], &(MAT1)->rows[row], &(MAT2)->rows[row], CTX); \
        dst->nnz += dst->rows[row].nnz; \
    } \
    return status; \
}

int gr_lil_mat_add(gr_lil_mat_t dst, const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)
{ GR_LIL_MAT_BOP(gr_sparse_vec_add, dst, mat1, mat2, ctx) }

int gr_lil_mat_sub(gr_lil_mat_t dst, const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)
{ GR_LIL_MAT_BOP(gr_sparse_vec_sub, dst, mat1, mat2, ctx) }

int gr_lil_mat_mul(gr_lil_mat_t dst, const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)
{ GR_LIL_MAT_BOP(gr_sparse_vec_mul, dst, mat1, mat2, ctx) }

#define GR_LIL_MAT_ACCUM_OP(FUNC, RES, MAT, C, CTX) \
{ \
    int row; \
    int status = GR_SUCCESS; \
                \
    if (gr_mat_is_compatible(RES, MAT, CTX) == T_FALSE) \
    { \
        return GR_DOMAIN; \
    } \
    (RES)->nnz = 0; \
    for (row = 0; row < (RES)->r; row++) \
    { \
        status |= FUNC(&(RES)->rows[row], &(MAT)->rows[row], C, CTX); \
        (RES)->nnz += (RES)->rows[row].nnz; \
    } \
    return status; \
}

int gr_lil_mat_addmul_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{
    GR_LIL_MAT_ACCUM_OP(gr_sparse_vec_addmul_scalar, dst, src, c, ctx)
}

int gr_lil_mat_submul_scalar(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{
    GR_LIL_MAT_ACCUM_OP(gr_sparse_vec_submul_scalar, dst, src, c, ctx)
}
