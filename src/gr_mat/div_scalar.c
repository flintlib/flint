/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_vec.h"

/* preinvert when appropriate */

#define GR_MAT_DIV_SCALAR(SCALAR_TYPE, res, mat, x, ctx) \
    slong i, r, c;                                       \
    int status = GR_SUCCESS;                             \
    r = gr_mat_nrows(res, ctx);                          \
    c = gr_mat_ncols(res, ctx);                          \
    if (c != 0)                                          \
        for (i = 0; i < r; i++)                          \
            status |= _gr_vec_div_##SCALAR_TYPE(         \
                res->rows[i],                            \
                mat->rows[i],                            \
                c, x, ctx                                \
            );                                           \
    return status;                                       \


int gr_mat_div_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
{ GR_MAT_DIV_SCALAR(scalar, res, mat, x, ctx); }

int gr_mat_div_scalar_si(gr_mat_t res, const gr_mat_t mat, slong x, gr_ctx_t ctx)
{ GR_MAT_DIV_SCALAR(scalar_si, res, mat, x, ctx); }

int gr_mat_div_scalar_ui(gr_mat_t res, const gr_mat_t mat, ulong x, gr_ctx_t ctx)
{ GR_MAT_DIV_SCALAR(scalar_ui, res, mat, x, ctx); }

int gr_mat_div_scalar_fmpz(gr_mat_t res, const gr_mat_t mat, fmpz_t x, gr_ctx_t ctx)
{ GR_MAT_DIV_SCALAR(scalar_fmpz, res, mat, x, ctx); }

int gr_mat_div_scalar_fmpq(gr_mat_t res, const gr_mat_t mat, fmpq_t x, gr_ctx_t ctx)
{ GR_MAT_DIV_SCALAR(scalar_fmpq, res, mat, x, ctx); }
