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


WARN_UNUSED_RESULT int gr_mat_update_lil_mat_nz(gr_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src->r; ++row)
        success |= gr_vec_update_sparse_vec_nz(dst->rows[row], &src->rows[row], ctx);
    return success;
}

WARN_UNUSED_RESULT int gr_mat_add_lil_mat(gr_mat_t dst, const gr_mat_t src1, const gr_lil_mat_t src2, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src1, ctx) != T_TRUE || gr_mat_is_compatible(dst, src2, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src1->r; ++row)
        success |= gr_vec_add_sparse_vec(dst->rows[row], src1->rows[row], &src2->rows[row], ctx);
    return success;
}

WARN_UNUSED_RESULT int gr_mat_sub_lil_mat(gr_mat_t dst, const gr_mat_t src1, const gr_lil_mat_t src2, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src1, ctx) != T_TRUE || gr_mat_is_compatible(dst, src2, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src1->r; ++row)
        success |= gr_vec_sub_sparse_vec(dst->rows[row], src1->rows[row], &src2->rows[row], ctx);
    return success;
}
WARN_UNUSED_RESULT int gr_mat_mul_lil_mat_nz(gr_mat_t dst, const gr_mat_t src1, const gr_lil_mat_t src2, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src1, ctx) != T_TRUE || gr_mat_is_compatible(dst, src2, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src1->r; ++row)
        success |= gr_vec_mul_sparse_vec_nz(dst->rows[row], src1->rows[row], &src2->rows[row], ctx);
    return success;
}
WARN_UNUSED_RESULT int gr_mat_div_lil_mat_nz(gr_mat_t dst, const gr_mat_t src1, const gr_lil_mat_t src2, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src1, ctx) != T_TRUE || gr_mat_is_compatible(dst, src2, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src1->r; ++row)
        success |= gr_vec_div_sparse_vec_nz(dst->rows[row], src1->rows[row], &src2->rows[row], ctx);
    return success;
}

WARN_UNUSED_RESULT int gr_mat_addmul_lil_mat_scalar(gr_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src->r; ++row)
        success |= gr_vec_addmul_sparse_vec_scalar(dst->rows[row], &src->rows[row], c, ctx);
    return success;
}

WARN_UNUSED_RESULT int gr_mat_submul_lil_mat_scalar(gr_mat_t dst, const gr_lil_mat_t src, gr_srcptr c, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src->r; ++row)
        success |= gr_vec_submul_sparse_vec_scalar(dst->rows[row], &src->rows[row], c, ctx);
    return success;
}

WARN_UNUSED_RESULT int gr_mat_addmul_lil_mat_scalar_si(gr_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src->r; ++row)
        success |= gr_vec_addmul_sparse_vec_scalar_si(dst->rows[row], &src->rows[row], c, ctx);
    return success;
}

WARN_UNUSED_RESULT int gr_mat_submul_lil_mat_scalar_si(gr_mat_t dst, const gr_lil_mat_t src, slong c, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src->r; ++row)
        success |= gr_vec_submul_sparse_vec_scalar_si(dst->rows[row], &src->rows[row], c, ctx);
    return success;
}

WARN_UNUSED_RESULT int gr_mat_addmul_lil_mat_scalar_fmpz(gr_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src->r; ++row)
        success |= gr_vec_addmul_sparse_vec_scalar_fmpz(dst->rows[row], &src->rows[row], c, ctx);
    return success;
}

WARN_UNUSED_RESULT int gr_mat_submul_lil_mat_scalar_fmpz(gr_mat_t dst, const gr_lil_mat_t src, const fmpz_t c, gr_ctx_t ctx)
{
    slong row;
    int success = GR_SUCCESS;
    if (gr_mat_is_compatible(dst, src, ctx) != T_TRUE)
        return GR_DOMAIN;
    for (row = 0; row < src->r; ++row)
        success |= gr_vec_submul_sparse_vec_scalar_fmpz(dst->rows[row], &src->rows[row], c, ctx);
    return success;
}
