/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_mat_lu_classical(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)
{
    acb_t d, e;
    acb_ptr * a;
    slong i, j, m, n, r, row, col;
    int result;

    if (acb_mat_is_empty(A))
        return 1;

    m = acb_mat_nrows(A);
    n = acb_mat_ncols(A);

    acb_mat_set(LU, A);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    acb_init(d);
    acb_init(e);

    result = 1;

    while (row < m && col < n)
    {
        r = acb_mat_find_pivot_partial(LU, row, m, col);

        if (r == -1)
        {
            result = 0;
            break;
        }
        else if (r != row)
            acb_mat_swap_rows(LU, P, row, r);

        acb_set(d, a[row] + col);

        for (j = row + 1; j < m; j++)
        {
            acb_div(e, a[j] + col, d, prec);
            acb_neg(e, e);
            _acb_vec_scalar_addmul(a[j] + col,
                a[row] + col, n - col, e, prec);
            acb_zero(a[j] + col);
            acb_neg(a[j] + row, e);
        }

        row++;
        col++;
    }

    acb_clear(d);
    acb_clear(e);

    return result;
}

