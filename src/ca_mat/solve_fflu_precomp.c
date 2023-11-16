/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_solve_fflu_precomp(ca_mat_t X, const slong * perm,
    const ca_mat_t A, const ca_t den, const ca_mat_t B, ca_ctx_t ctx)
{
    ca_t t;
    slong i, j, k, c, n, m;

    n = ca_mat_nrows(X);
    m = ca_mat_ncols(X);

    if (X == B)
    {
        ca_ptr tmp = flint_malloc(sizeof(ca_struct) * n);

        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
                tmp[i] = B->rows[perm[i]][c];
            for (i = 0; i < n; i++)
                X->rows[i][c] = tmp[i];
        }

        flint_free(tmp);
    }
    else
    {
        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
            {
                ca_set(ca_mat_entry(X, i, c),
                    ca_mat_entry(B, perm[i], c), ctx);
            }
        }
    }

    ca_init(t, ctx);

    for (k = 0; k < m; k++)
    {
        /* Fraction-free forward substitution */
        for (i = 0; i < n - 1; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                ca_mul(ca_mat_entry(X, j, k), ca_mat_entry(X, j, k), ca_mat_entry(A, i, i), ctx);
                ca_mul(t, ca_mat_entry(A, j, i), ca_mat_entry(X, i, k), ctx);
                ca_sub(ca_mat_entry(X, j, k), ca_mat_entry(X, j, k), t, ctx);
                if (i > 0)
                    ca_div(ca_mat_entry(X, j, k), ca_mat_entry(X, j, k), ca_mat_entry(A, i-1, i-1), ctx);
            }
        }

        /* Fraction-free back substitution */
        for (i = n - 2; i >= 0; i--)
        {
            ca_mul(ca_mat_entry(X, i, k), ca_mat_entry(X, i, k), ca_mat_entry(A, n-1, n-1), ctx);
            for (j = i + 1; j < n; j++)
            {
                ca_mul(t, ca_mat_entry(X, j, k), ca_mat_entry(A, i, j), ctx);
                ca_sub(ca_mat_entry(X, i, k), ca_mat_entry(X, i, k), t, ctx);
            }
            ca_div(ca_mat_entry(X, i, k), ca_mat_entry(X, i, k), ca_mat_entry(A, i, i), ctx);
        }
    }

    ca_mat_div_ca(X, X, den, ctx);

    ca_clear(t, ctx);
}
