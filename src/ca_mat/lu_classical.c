/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

int
ca_mat_lu_classical(slong * res_rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
{
    ca_t d, e;
    ca_ptr * a;
    slong i, j, m, n, r, rank, row, col;
    int success;
    truth_t found_pivot;

    if (ca_mat_is_empty(A))
    {
        *res_rank = 0;
        return 1;
    }

    m = ca_mat_nrows(A);
    n = ca_mat_ncols(A);

    ca_mat_set(LU, A, ctx);

    a = LU->rows;

    rank = row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    ca_init(d, ctx);
    ca_init(e, ctx);

    success = 1;

    while (row < m && col < n)
    {
        found_pivot = ca_mat_find_pivot(&r, LU, row, m, col, ctx);

        if (found_pivot == T_UNKNOWN)
        {
            success = 0;
            break;
        }

        if (found_pivot == T_FALSE)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }
            col++;
            continue;
        }

        rank++;

        if (r != row)
            _ca_mat_swap_rows(LU, P, row, r);

        ca_inv(d, a[row] + col, ctx);

        for (j = row + 1; j < m; j++)
        {
            ca_mul(e, a[j] + col, d, ctx);
            ca_neg(e, e, ctx);

            _ca_vec_scalar_addmul_ca(a[j] + col + 1, a[row] + col + 1, n - col - 1, e, ctx);

            ca_zero(a[j] + col, ctx);
            ca_neg(a[j] + rank - 1, e, ctx);
        }

        row++;
        col++;
    }

    ca_clear(d, ctx);
    ca_clear(e, ctx);

    *res_rank = rank;
    return success;
}
