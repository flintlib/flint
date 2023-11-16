/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
arb_mat_lu_classical(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec)
{
    arb_t d, e;
    arb_ptr * a;
    slong i, j, m, n, r, row, col;
    int result;

    if (arb_mat_is_empty(A))
        return 1;

    m = arb_mat_nrows(A);
    n = arb_mat_ncols(A);

    arb_mat_set(LU, A);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    arb_init(d);
    arb_init(e);

    result = 1;

    while (row < m && col < n)
    {
        r = arb_mat_find_pivot_partial(LU, row, m, col);

        if (r == -1)
        {
            result = 0;
            break;
        }
        else if (r != row)
            arb_mat_swap_rows(LU, P, row, r);

        arb_set(d, a[row] + col);

        for (j = row + 1; j < m; j++)
        {
            arb_div(e, a[j] + col, d, prec);
            arb_neg(e, e);
            _arb_vec_scalar_addmul(a[j] + col,
                a[row] + col, n - col, e, prec);
            arb_zero(a[j] + col);
            arb_neg(a[j] + row, e);
        }

        row++;
        col++;
    }

    arb_clear(d);
    arb_clear(e);

    return result;
}

