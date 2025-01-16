/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"


static inline int
TEMPLATE(T, mat_pivot) (TEMPLATE(T, mat_t) A, slong * P, slong start_row,
                        slong col, const TEMPLATE(T, ctx_t) ctx)
{
    slong j;

    if (!TEMPLATE(T, is_zero)
        (TEMPLATE(T, mat_entry) (A, start_row, col), ctx))
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (!TEMPLATE(T, is_zero) (TEMPLATE(T, mat_entry) (A, j, col), ctx))
        {
            TEMPLATE(T, mat_swap_rows) (A, P, j, start_row, ctx);
            return -1;
        }
    }
    return 0;
}


slong
TEMPLATE(T, mat_lu_classical) (slong * P,
                               TEMPLATE(T, mat_t) A,
                               int rank_check, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) d, e, neg_e;
    slong i, m, n, rank, length, row, col;

    m = A->r;
    n = A->c;

    rank = row = col = 0;

    for (i = 0; i < m; i++)
        P[i] = i;

    TEMPLATE(T, init) (d, ctx);
    TEMPLATE(T, init) (e, ctx);
    TEMPLATE(T, init) (neg_e, ctx);

    while (row < m && col < n)
    {
        if (TEMPLATE(T, mat_pivot) (A, P, row, col, ctx) == 0)
        {
            if (rank_check)
            {
                rank = 0;
                goto cleanup;
            }
            col++;
            continue;
        }

        rank++;

        TEMPLATE(T, inv) (d, TEMPLATE(T, mat_entry) (A, row, col), ctx);

        length = n - col - 1;

        for (i = row + 1; i < m; i++)
        {
            TEMPLATE(T, mul) (e, TEMPLATE(T, mat_entry) (A, i, col), d, ctx);
            if (length != 0)
            {
                TEMPLATE(T, neg) (neg_e, e, ctx);
                _TEMPLATE3(T, vec_scalar_addmul, T) (TEMPLATE(T, mat_entry) (A, i, col + 1),
                                                     TEMPLATE(T, mat_entry) (A, row, col + 1),
                                                     length, neg_e, ctx);
            }

            TEMPLATE(T, zero) (TEMPLATE(T, mat_entry) (A, i, col), ctx);
            TEMPLATE(T, set) (TEMPLATE(T, mat_entry) (A, i, rank - 1), e, ctx);
        }
        row++;
        col++;
    }

cleanup:
    TEMPLATE(T, clear) (d, ctx);
    TEMPLATE(T, clear) (e, ctx);
    TEMPLATE(T, clear) (neg_e, ctx);

    return rank;
}


#endif
