/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"

static inline int
fmpz_mod_mat_pivot(fmpz_mod_mat_t A, slong * P, slong start_row,
                        slong col, const fmpz_mod_ctx_t ctx)
{
    slong j, t;
    fmpz * u;

    if (!fmpz_is_zero(fmpz_mod_mat_entry(A, start_row, col)))
        return 1;

    for (j = start_row + 1; j < A->mat->r; j++)
    {
        if (!fmpz_is_zero(fmpz_mod_mat_entry(A, j, col)))
        {
            u = A->mat->rows[j];
            A->mat->rows[j] = A->mat->rows[start_row];
            A->mat->rows[start_row] = u;

            t = P[j];
            P[j] = P[start_row];
            P[start_row] = t;

            return -1;
        }
    }
    return 0;
}


slong fmpz_mod_mat_lu_classical(slong * P, fmpz_mod_mat_t A, int rank_check)
{
    fmpz_t d, e, neg_e;
    fmpz ** a;
    slong i, m, n, rank, length, row, col;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init(ctx, A->mod);

    m = A->mat->r;
    n = A->mat->c;
    a = A->mat->rows;

    rank = row = col = 0;

    for (i = 0; i < m; i++)
        P[i] = i;

    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(neg_e);

    while (row < m && col < n)
    {
        if (fmpz_mod_mat_pivot(A, P, row, col, ctx) == 0)
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

        fmpz_mod_inv(d, a[row] + col, ctx);

        length = n - col - 1;

        for (i = row + 1; i < m; i++)
        {
            fmpz_mod_mul(e, a[i] + col, d, ctx);
            if (length != 0)
            {
                fmpz_mod_neg(neg_e, e, ctx);
                _fmpz_mod_vec_scalar_addmul_fmpz_mod(a[i] + col + 1,
                                                     a[row] + col + 1,
                                                     length, neg_e, ctx);
            }

            fmpz_zero(a[i] + col);
            fmpz_set(a[i] + rank - 1, e);
        }
        row++;
        col++;
    }

cleanup:
    fmpz_clear(d);
    fmpz_clear(e);
    fmpz_clear(neg_e);

    fmpz_mod_ctx_clear(ctx);

    return rank;
}

