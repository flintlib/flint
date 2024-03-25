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

int gr_csr_mat_neg(gr_csr_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
    {
        return GR_DOMAIN;
    }
    status |= gr_csr_mat_set(dst, src, ctx);
    status |= _gr_vec_neg(dst->nzs, dst->nzs, dst->nnz, ctx);
    return status;
}

int gr_lil_mat_neg(gr_lil_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    int row;
    int status = GR_SUCCESS;

    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
    {
        return GR_DOMAIN;
    }
    dst->nnz = src->nnz;
    for (row = 0; row < dst->r; row++)
{
        status |= gr_sparse_vec_neg(&dst->rows[row], &src->rows[row], ctx);
    }
    return status;
}

int gr_coo_mat_neg(gr_coo_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (gr_mat_is_compatible(dst, src, ctx) == T_FALSE)
    {
        return GR_DOMAIN;
    }
    status |= gr_coo_mat_set(dst, src, ctx);
    status |= _gr_vec_neg(dst->nzs, dst->nzs, dst->nnz, ctx);
    return status;
}
