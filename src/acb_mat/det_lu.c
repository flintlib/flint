/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

slong
acb_mat_gauss_partial(acb_mat_t A, slong prec)
{
    acb_t e;
    acb_ptr * a;
    slong j, m, n, r, rank, row, col, sign;

    m = A->r;
    n = A->c;
    a = A->rows;
    rank = row = col = 0;
    sign = 1;

    acb_init(e);

    while (row < m && col < n)
    {
        r = acb_mat_find_pivot_partial(A, row, m, col);

        if (r == -1)
        {
            break;
        }
        else if (r != row)
        {
            acb_mat_swap_rows(A, NULL, row, r);
            sign *= -1;
        }

        rank++;

        for (j = row + 1; j < m; j++)
        {
            acb_div(e, a[j] + col, a[row] + col, prec);
            acb_neg(e, e);
            _acb_vec_scalar_addmul(a[j] + col + 1, a[row] + col + 1, n - col - 1, e, prec);
        }

        row++;
        col++;
    }

    acb_clear(e);

    return rank * sign;
}

void
acb_vec_get_arf_2norm_squared_bound(arf_t s, acb_srcptr vec, slong len, slong prec)
{
    slong i;
    arf_t t;

    arf_init(t);
    arf_zero(s);

    for (i = 0; i < len; i++)
    {
        arb_get_abs_ubound_arf(t, acb_realref(vec + i), prec);
        arf_addmul(s, t, t, prec, ARF_RND_UP);
        arb_get_abs_ubound_arf(t, acb_imagref(vec + i), prec);
        arf_addmul(s, t, t, prec, ARF_RND_UP);
    }

    arf_clear(t);
}

void
acb_mat_det_lu_inplace(acb_t det, acb_mat_t A, slong prec)
{
    slong i, n, sign, rank;
    int is_real;

    n = acb_mat_nrows(A);
    rank = acb_mat_gauss_partial(A, prec);
    sign = (rank < 0) ? -1 : 1;
    rank = FLINT_ABS(rank);

    _acb_mat_diag_prod(det, A, 0, rank, prec);
    acb_mul_si(det, det, sign, prec);

    /* bound unreduced part using Hadamard's inequality */
    if (rank < n)
    {
        arf_t t;
        arf_t d;
        acb_t e;

        arf_init(t);
        arf_init(d);
        acb_init(e);

        arf_one(d);

        is_real = acb_mat_is_real(A);

        for (i = rank; i < n; i++)
        {
            acb_vec_get_arf_2norm_squared_bound(t, A->rows[i] + rank,  n - rank, MAG_BITS);
            arf_mul(d, d, t, MAG_BITS, ARF_RND_UP);
        }

        /* now d contains the absolute value of the determinant */
        arf_sqrt(d, d, MAG_BITS, ARF_RND_UP);

        /* multiply by disc with radius d */
        if (is_real)
        {
            arb_add_error_arf(acb_realref(e), d);
        }
        else
        {
            arb_add_error_arf(acb_realref(e), d);
            arb_add_error_arf(acb_imagref(e), d);
        }

        acb_mul(det, det, e, prec);

        acb_clear(e);
        arf_clear(d);
        arf_clear(t);
    }
}

void
acb_mat_det_lu(acb_t det, const acb_mat_t A, slong prec)
{
    slong n;

    n = acb_mat_nrows(A);

    if (n == 0)
    {
        acb_one(det);
    }
    else
    {
        acb_mat_t T;
        acb_mat_init(T, acb_mat_nrows(A), acb_mat_ncols(A));
        acb_mat_set(T, A);
        acb_mat_det_lu_inplace(det, T, prec);
        acb_mat_clear(T);
    }
}

