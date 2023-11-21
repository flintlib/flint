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
ca_mat_fflu(slong * res_rank, slong * P, ca_mat_t LU, ca_t den, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
{
    ca_t d, e;
    slong i, j, k, m, n, r, rank, row, col;
    int success;
    truth_t found_pivot;

    if (ca_mat_is_empty(A))
    {
        *res_rank = 0;
        ca_one(den, ctx);
        return 1;
    }

    m = ca_mat_nrows(A);
    n = ca_mat_ncols(A);

    ca_mat_set(LU, A, ctx);

    rank = row = col = 0;
    if (P != NULL)
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
                ca_zero(den, ctx);
                rank = 0;
                break;
            }
            col++;
            continue;
        }

        rank++;

        if (r != row)
            _ca_mat_swap_rows(LU, P, row, r);

        if (row > 0)
            ca_inv(d, den, ctx);

        for (j = row + 1; j < m; j++)
        {
            for (k = col + 1; k < n; k++)
            {
                ca_mul(ca_mat_entry(LU, j, k), ca_mat_entry(LU, j, k), ca_mat_entry(LU, row, col), ctx);
                ca_mul(e, ca_mat_entry(LU, j, col), ca_mat_entry(LU, row, k), ctx);
                ca_sub(ca_mat_entry(LU, j, k), ca_mat_entry(LU, j, k), e, ctx);
                if (row > 0)
                    ca_mul(ca_mat_entry(LU, j, k), ca_mat_entry(LU, j, k), d, ctx);
            }

            /* todo: zero at (j, col) ? */
        }

        ca_set(den, ca_mat_entry(LU, row, col), ctx);
        row++;
        col++;
    }

    ca_clear(d, ctx);
    ca_clear(e, ctx);

    if (success)
    {
        if (rank == 0)
            ca_zero(den, ctx);
    }

    if (!success)
        ca_unknown(den, ctx);

    *res_rank = rank;
    return success;
}
