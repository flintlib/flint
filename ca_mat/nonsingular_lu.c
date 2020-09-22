/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

truth_t
ca_mat_nonsingular_lu(slong * P, ca_mat_t LU, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_t d, e;
    ca_ptr * a;
    slong i, j, m, n, r, row, col;
    truth_t result;

    if (ca_mat_is_empty(A))
        return T_TRUE;

    m = arb_mat_nrows(A);
    n = arb_mat_ncols(A);

    ca_mat_set(LU, A, ctx);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    ca_init(d, ctx);
    ca_init(e, ctx);

    result = T_TRUE;

    while (row < m && col < n)
    {
        result = ca_mat_find_pivot(&r, LU, row, m, col, ctx);

        if (result != T_TRUE)
            break;

        if (r != row)
            _ca_mat_swap_rows(LU, P, row, r);

        ca_inv(d, a[row] + col, ctx);

        for (j = row + 1; j < m; j++)
        {
            ca_mul(e, a[j] + col, d, ctx);
            ca_neg(e, e, ctx);
            ca_vec_scalar_addmul_ca(a[j] + col, a[row] + col, n - col, e, ctx);
            ca_zero(a[j] + col, ctx);
            ca_neg(a[j] + row, e, ctx);
        }

        row++;
        col++;
    }

    ca_clear(d, ctx);
    ca_clear(e, ctx);

    return result;
}
