/*
    Copyright (C) 2018 arbguest
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_approx_solve_lu_precomp(arb_mat_t X, const slong * perm,
    const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong i, c, n, m;

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

    arb_mat_get_mid(X, X);
    arb_mat_approx_solve_tril(X, A, X, 1, prec);
    arb_mat_approx_solve_triu(X, A, X, 0, prec);
}
