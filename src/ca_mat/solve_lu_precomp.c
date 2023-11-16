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
ca_mat_solve_lu_precomp(ca_mat_t X, const slong * perm,
    const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
{
    slong i, c, n, m;

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

    /* todo: inline for small n? */
    ca_mat_solve_tril(X, A, X, 1, ctx);
    ca_mat_solve_triu(X, A, X, 0, ctx);
}
