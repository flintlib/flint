/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

static void _stirling_number_1u_vec_next(ca_ptr row,
    ca_srcptr prev, slong n, slong klen, ca_ctx_t ctx)
{
    slong k;

    if (klen > n) ca_one(row + n, ctx);
    if (n != 0 && klen != 0) ca_zero(row, ctx);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        ca_mul_ui(row + k, prev + k, n - 1, ctx);
        ca_add(row + k, prev + k - 1, row + k, ctx);
    }

    for (k = n + 1; k < klen; k++)
        ca_zero(row + k, ctx);
}

static void _stirling_number_1_vec_next(ca_ptr row,
    ca_srcptr prev, slong n, slong klen, ca_ctx_t ctx)
{
    slong k;

    if (klen > n) ca_one(row + n, ctx);
    if (n != 0 && klen != 0) ca_zero(row, ctx);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        ca_mul_ui(row + k, prev + k, n - 1, ctx);
        ca_sub(row + k, prev + k - 1, row + k, ctx);
    }

    for (k = n + 1; k < klen; k++)
        ca_zero(row + k, ctx);
}

static void _stirling_number_2_vec_next(ca_ptr row,
    ca_srcptr prev, slong n, slong klen, ca_ctx_t ctx)
{
    slong k;

    if (klen > n) ca_one(row + n, ctx);
    if (n != 0 && klen != 0) ca_zero(row, ctx);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        ca_mul_ui(row + k, prev + k, k, ctx);
        ca_add(row + k, prev + k - 1, row + k, ctx);
    }

    for (k = n + 1; k < klen; k++)
        ca_zero(row + k, ctx);
}

static void
_stirling_matrix_1u(ca_mat_t mat, ca_ctx_t ctx)
{
    slong n;

    if (ca_mat_is_empty(mat))
        return;

    for (n = 0; n < mat->r; n++)
        _stirling_number_1u_vec_next(mat->rows[n],
            mat->rows[n - (n != 0)], n, mat->c, ctx);
}

static void
_stirling_matrix_1(ca_mat_t mat, ca_ctx_t ctx)
{
    slong n;

    if (ca_mat_is_empty(mat))
        return;

    for (n = 0; n < mat->r; n++)
        _stirling_number_1_vec_next(mat->rows[n],
            mat->rows[n - (n != 0)], n, mat->c, ctx);
}

static void
_stirling_matrix_2(ca_mat_t mat, ca_ctx_t ctx)
{
    slong n;

    if (ca_mat_is_empty(mat))
        return;

    for (n = 0; n < mat->r; n++)
        _stirling_number_2_vec_next(mat->rows[n],
            mat->rows[n - (n != 0)], n, mat->c, ctx);
}

void
ca_mat_stirling(ca_mat_t mat, int kind, ca_ctx_t ctx)
{
    if (kind == 0)
        _stirling_matrix_1u(mat, ctx);
    else if (kind == 1)
        _stirling_matrix_1(mat, ctx);
    else
        _stirling_matrix_2(mat, ctx);
}
