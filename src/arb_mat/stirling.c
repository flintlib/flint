/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

static void _stirling_number_1u_vec_next(arb_ptr row,
    arb_srcptr prev, slong n, slong klen, slong prec)
{
    slong k;

    if (klen > n) arb_one(row + n);
    if (n != 0 && klen != 0) arb_zero(row);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        arb_mul_ui(row + k, prev + k, n - 1, prec);
        arb_add(row + k, prev + k - 1, row + k, prec);
    }

    for (k = n + 1; k < klen; k++)
        arb_zero(row + k);
}

static void _stirling_number_1_vec_next(arb_ptr row,
    arb_srcptr prev, slong n, slong klen, slong prec)
{
    slong k;

    if (klen > n) arb_one(row + n);
    if (n != 0 && klen != 0) arb_zero(row);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        arb_mul_ui(row + k, prev + k, n - 1, prec);
        arb_sub(row + k, prev + k - 1, row + k, prec);
    }

    for (k = n + 1; k < klen; k++)
        arb_zero(row + k);
}

static void _stirling_number_2_vec_next(arb_ptr row,
    arb_srcptr prev, slong n, slong klen, slong prec)
{
    slong k;

    if (klen > n) arb_one(row + n);
    if (n != 0 && klen != 0) arb_zero(row);

    for (k = FLINT_MIN(n, klen) - 1; k >= 1; k--)
    {
        arb_mul_ui(row + k, prev + k, k, prec);
        arb_add(row + k, prev + k - 1, row + k, prec);
    }

    for (k = n + 1; k < klen; k++)
        arb_zero(row + k);
}

static void
_stirling_matrix_1u(arb_mat_t mat, slong prec)
{
    slong n;

    if (arb_mat_is_empty(mat))
        return;

    for (n = 0; n < mat->r; n++)
        _stirling_number_1u_vec_next(mat->rows[n],
            mat->rows[n - (n != 0)], n, mat->c, prec);
}

static void
_stirling_matrix_1(arb_mat_t mat, slong prec)
{
    slong n;

    if (arb_mat_is_empty(mat))
        return;

    for (n = 0; n < mat->r; n++)
        _stirling_number_1_vec_next(mat->rows[n],
            mat->rows[n - (n != 0)], n, mat->c, prec);
}

static void
_stirling_matrix_2(arb_mat_t mat, slong prec)
{
    slong n;

    if (arb_mat_is_empty(mat))
        return;

    for (n = 0; n < mat->r; n++)
        _stirling_number_2_vec_next(mat->rows[n],
            mat->rows[n - (n != 0)], n, mat->c, prec);
}

void
arb_mat_stirling(arb_mat_t mat, int kind, slong prec)
{
    if (kind == 0)
        _stirling_matrix_1u(mat, prec);
    else if (kind == 1)
        _stirling_matrix_1(mat, prec);
    else
        _stirling_matrix_2(mat, prec);
}

