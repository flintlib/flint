/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

slong
arb_mat_gauss_partial(arb_mat_t A, slong prec)
{
    arb_t e;
    arb_ptr * a;
    slong j, m, n, r, rank, row, col, sign;

    m = A->r;
    n = A->c;
    a = A->rows;
    rank = row = col = 0;
    sign = 1;

    arb_init(e);

    while (row < m && col < n)
    {
        r = arb_mat_find_pivot_partial(A, row, m, col);

        if (r == -1)
        {
            break;
        }
        else if (r != row)
        {
            arb_mat_swap_rows(A, NULL, row, r);
            sign *= -1;
        }

        rank++;

        for (j = row + 1; j < m; j++)
        {
            arb_div(e, a[j] + col, a[row] + col, prec);
            arb_neg(e, e);
            _arb_vec_scalar_addmul(a[j] + col + 1, a[row] + col + 1, n - col - 1, e, prec);
        }

        row++;
        col++;
    }

    arb_clear(e);

    return rank * sign;
}

void
arb_vec_get_arf_2norm_squared_bound(arf_t s, arb_srcptr vec, slong len, slong prec)
{
    slong i;
    arf_t t;

    arf_init(t);
    arf_zero(s);

    for (i = 0; i < len; i++)
    {
        arb_get_abs_ubound_arf(t, vec + i, prec);
        arf_addmul(s, t, t, prec, ARF_RND_UP);
    }

    arf_clear(t);
}

void
arb_mat_det_lu_inplace(arb_t det, arb_mat_t A, slong prec)
{
    slong i, n, sign, rank;

    n = arb_mat_nrows(A);
    rank = arb_mat_gauss_partial(A, prec);
    sign = (rank < 0) ? -1 : 1;
    rank = FLINT_ABS(rank);

    arb_set_si(det, sign);
    for (i = 0; i < rank; i++)
        arb_mul(det, det, arb_mat_entry(A, i, i), prec);

    /* bound unreduced part using Hadamard's inequality */
    if (rank < n)
    {
        arf_t t;
        arf_t d;
        arb_t b;

        arf_init(t);
        arf_init(d);
        arb_init(b);

        arf_one(d);

        for (i = rank; i < n; i++)
        {
            arb_vec_get_arf_2norm_squared_bound(t, A->rows[i] + rank,  n - rank, MAG_BITS);
            arf_mul(d, d, t, MAG_BITS, ARF_RND_UP);
        }

        arf_sqrt(d, d, MAG_BITS, ARF_RND_UP);
        arb_add_error_arf(b, d);

        arb_mul(det, det, b, prec);

        arf_clear(d);
        arf_clear(t);
        arb_clear(b);
    }
}

void
arb_mat_det_lu(arb_t det, const arb_mat_t A, slong prec)
{
    slong n;

    n = arb_mat_nrows(A);

    if (n == 0)
    {
        arb_one(det);
    }
    else
    {
        arb_mat_t T;
        arb_mat_init(T, arb_mat_nrows(A), arb_mat_ncols(A));
        arb_mat_set(T, A);
        arb_mat_det_lu_inplace(det, T, prec);
        arb_mat_clear(T);
    }
}
