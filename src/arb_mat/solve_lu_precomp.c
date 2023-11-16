/*
    Copyright (C) 2012,2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_solve_lu_precomp(arb_mat_t X, const slong * perm,
    const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong i, j, c, n, m;

    n = arb_mat_nrows(X);
    m = arb_mat_ncols(X);

    if (X == B)
    {
        arb_ptr tmp = flint_malloc(sizeof(arb_struct) * n);

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
                arb_set(arb_mat_entry(X, i, c),
                    arb_mat_entry(B, perm[i], c));
            }
        }
    }

    /* solve_tril and solve_triu have some overhead */
    if (n >= 4)
    {
        arb_mat_solve_tril(X, A, X, 1, prec);
        arb_mat_solve_triu(X, A, X, 0, prec);
        return;
    }

    for (c = 0; c < m; c++)
    {
        /* solve Ly = b */
        for (i = 1; i < n; i++)
        {
            for (j = 0; j < i; j++)
            {
                arb_submul(arb_mat_entry(X, i, c),
                    arb_mat_entry(A, i, j), arb_mat_entry(X, j, c), prec);
            }
        }

        /* solve Ux = y */
        for (i = n - 1; i >= 0; i--)
        {
            for (j = i + 1; j < n; j++)
            {
                arb_submul(arb_mat_entry(X, i, c),
                    arb_mat_entry(A, i, j), arb_mat_entry(X, j, c), prec);
            }

            arb_div(arb_mat_entry(X, i, c), arb_mat_entry(X, i, c),
                arb_mat_entry(A, i, i), prec);
        }
    }
}
