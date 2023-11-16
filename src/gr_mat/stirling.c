/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/* todo: optimize for fmpz */

#define ENTRY(vec, i) GR_ENTRY(vec, i, sz)

static int _stirling_number_1u_vec_next(gr_ptr row,
    gr_srcptr prev, slong n, slong klen, gr_ctx_t ctx)
{
    slong k;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (klen > n)
        status |= gr_one(ENTRY(row, n), ctx);

    if (n != 0 && klen != 0)
        status |= gr_zero(row, ctx);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        status |= gr_mul_ui(ENTRY(row, k), ENTRY(prev, k), n - 1, ctx);
        status |= gr_add(ENTRY(row, k), ENTRY(prev, k - 1), ENTRY(row, k), ctx);
    }

    for (k = n + 1; k < klen; k++)
        status |= gr_zero(ENTRY(row, k), ctx);

    return status;
}

static int _stirling_number_1_vec_next(gr_ptr row,
    gr_srcptr prev, slong n, slong klen, gr_ctx_t ctx)
{
    slong k;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (klen > n)
        status |= gr_one(ENTRY(row, n), ctx);

    if (n != 0 && klen != 0)
        status |= gr_zero(row, ctx);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        status |= gr_mul_ui(ENTRY(row, k), ENTRY(prev, k), n - 1, ctx);
        status |= gr_sub(ENTRY(row, k), ENTRY(prev, k - 1), ENTRY(row, k), ctx);
    }

    for (k = n + 1; k < klen; k++)
        status |= gr_zero(ENTRY(row, k), ctx);

    return status;
}

static int _stirling_number_2_vec_next(gr_ptr row,
    gr_srcptr prev, slong n, slong klen, gr_ctx_t ctx)
{
    slong k;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (klen > n)
        status |= gr_one(ENTRY(row, n), ctx);

    if (n != 0 && klen != 0)
        status |= gr_zero(row, ctx);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        status |= gr_mul_ui(ENTRY(row, k), ENTRY(prev, k), k, ctx);
        status |= gr_add(ENTRY(row, k), ENTRY(prev, k - 1), ENTRY(row, k), ctx);
    }

    for (k = n + 1; k < klen; k++)
        status |= gr_zero(ENTRY(row, k), ctx);

    return status;
}

static int
_stirling_matrix_1u(gr_mat_t mat, gr_ctx_t ctx)
{
    slong n;
    int status = GR_SUCCESS;

    if (mat->c != 0)
        for (n = 0; n < mat->r; n++)
            status |= _stirling_number_1u_vec_next(mat->rows[n],
                mat->rows[n - (n != 0)], n, mat->c, ctx);

    return status;
}

static int
_stirling_matrix_1(gr_mat_t mat, gr_ctx_t ctx)
{
    slong n;
    int status = GR_SUCCESS;

    if (mat->c != 0)
        for (n = 0; n < mat->r; n++)
            status |= _stirling_number_1_vec_next(mat->rows[n],
                mat->rows[n - (n != 0)], n, mat->c, ctx);

    return status;
}

static int
_stirling_matrix_2(gr_mat_t mat, gr_ctx_t ctx)
{
    slong n;
    int status = GR_SUCCESS;

    if (mat->c != 0)
        for (n = 0; n < mat->r; n++)
            status |= _stirling_number_2_vec_next(mat->rows[n],
                mat->rows[n - (n != 0)], n, mat->c, ctx);

    return status;
}

int
gr_mat_stirling(gr_mat_t mat, int kind, gr_ctx_t ctx)
{
    if (kind == 0)
        return _stirling_matrix_1u(mat, ctx);
    else if (kind == 1)
        return _stirling_matrix_1(mat, ctx);
    else if (kind == 2)
        return _stirling_matrix_2(mat, ctx);
    else
        return GR_DOMAIN;
}
