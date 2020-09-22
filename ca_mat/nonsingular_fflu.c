/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

#define E(i,j) ca_mat_entry(LU, i, j)

truth_t
ca_mat_nonsingular_fflu(slong * P, ca_mat_t LU, ca_t den, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_t e, t;
    ca_ptr * a;
    slong i, j, k, m, n, r, row, col;
    truth_t result;

    if (ca_mat_is_empty(A))
    {
        ca_one(den, ctx);
        return T_TRUE;
    }

    m = arb_mat_nrows(A);
    n = arb_mat_ncols(A);

    ca_mat_set(LU, A, ctx);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    ca_init(e, ctx);
    ca_init(t, ctx);

    result = T_TRUE;

    ca_one(den, ctx);

    while (row < m && col < n)
    {
        result = ca_mat_find_pivot(&r, LU, row, m, col, ctx);

        if (result != T_TRUE)
            break;

        if (r != row)
            _ca_mat_swap_rows(LU, P, row, r);

        ca_inv(e, den, ctx);

        for (j = row + 1; j < m; j++)
        {
            for (k = col + 1; k < n; k++)
            {
                ca_mul(E(j, k), E(j, k), E(row, col), ctx);
                ca_mul(t, E(j, col), E(row, k), ctx);
                ca_sub(E(j, k), E(j, k), t, ctx);
                if (row > 0)
                    ca_mul(E(j, k), E(j, k), e, ctx);
            }
        }

        ca_set(den, E(row, col), ctx);

        row++;
        col++;
    }

    if (result == T_FALSE)
        ca_zero(den, ctx);
    else if (result == T_UNKNOWN)
        ca_unknown(den, ctx);

    ca_clear(e, ctx);
    ca_clear(t, ctx);

    return result;
}

#undef E
