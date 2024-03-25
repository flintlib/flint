/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

#define GR_MAT_ADD_SCALAR(FUNC, res, mat, x, ctx)        \
    slong i, j, r, c, sz = (ctx)->sizeof_elem;           \
    int status = GR_SUCCESS;                             \
    r = gr_mat_nrows(res, ctx);                          \
    c = gr_mat_ncols(res, ctx);                          \
    if (res == mat)                                      \
    {                                                    \
        for (i = 0; i < FLINT_MIN(r, c); i++)            \
            status |= (FUNC)(                            \
                GR_MAT_ENTRY(res, i, i, sz),             \
                GR_MAT_ENTRY(res, i, i, sz),             \
                x, ctx                                   \
            );                                           \
    }                                                    \
    else                                                 \
    {                                                    \
        for (i = 0; i < r; i++)                          \
        {                                                \
            for (j = 0; j < c; j++)                      \
            {                                            \
                /* todo: vectorize */                    \
                if (i == j)                              \
                    status |= (FUNC)(                    \
                        GR_MAT_ENTRY(res, i, j, sz),     \
                        GR_MAT_ENTRY(mat, i, j, sz),     \
                        x, ctx                           \
                    );                                   \
                else                                     \
                    status |= gr_set(                    \
                        GR_MAT_ENTRY(res, i, j, sz),     \
                        GR_MAT_ENTRY(mat, i, j, sz),     \
                        ctx                              \
                    );                                   \
            }                                            \
        }                                                \
    }                                                    \
    return status;                                       \


int gr_mat_add_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
{ GR_MAT_ADD_SCALAR(gr_add, res, mat, x, ctx) }

int gr_mat_add_scalar_si(gr_mat_t res, const gr_mat_t mat, slong x, gr_ctx_t ctx)
{ GR_MAT_ADD_SCALAR(gr_add_si, res, mat, x, ctx) }

int gr_mat_add_scalar_ui(gr_mat_t res, const gr_mat_t mat, ulong x, gr_ctx_t ctx)
{ GR_MAT_ADD_SCALAR(gr_add_ui, res, mat, x, ctx) }

int gr_mat_add_scalar_fmpz(gr_mat_t res, const gr_mat_t mat, fmpz_t x, gr_ctx_t ctx)
{ GR_MAT_ADD_SCALAR(gr_add_fmpz, res, mat, x, ctx) }

int gr_mat_add_scalar_fmpq(gr_mat_t res, const gr_mat_t mat, fmpq_t x, gr_ctx_t ctx)
{ GR_MAT_ADD_SCALAR(gr_add_fmpq, res, mat, x, ctx) }
