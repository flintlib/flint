/*
    Copyright (C) 2018 arbguest

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_approx_solve_lu_precomp(acb_mat_t X, const slong * perm,
    const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong i, c, n, m;

    n = acb_mat_nrows(X);
    m = acb_mat_ncols(X);

    if (X == B)
    {
        acb_ptr tmp = flint_malloc(sizeof(acb_struct) * n);

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
                acb_set(acb_mat_entry(X, i, c),
                    acb_mat_entry(B, perm[i], c));
            }
        }
    }

    acb_mat_get_mid(X, X);
    acb_mat_approx_solve_tril(X, A, X, 1, prec);
    acb_mat_approx_solve_triu(X, A, X, 0, prec);
}
