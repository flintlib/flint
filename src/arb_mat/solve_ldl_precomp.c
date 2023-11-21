/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_solve_ldl_precomp(arb_mat_t X,
    const arb_mat_t L, const arb_mat_t B, slong prec)
{
    slong i, j, c, n, m;

    n = arb_mat_nrows(X);
    m = arb_mat_ncols(X);

    arb_mat_set(X, B);

    for (c = 0; c < m; c++)
    {
        /* solve Lz = b */
        for (i = 1; i < n; i++)
        {
            for (j = 0; j < i; j++)
            {
                arb_submul(arb_mat_entry(X, i, c),
                    arb_mat_entry(L, i, j), arb_mat_entry(X, j, c), prec);
            }
        }

        /* solve Dy = z */
        for (i = 0; i < n; i++)
        {
            arb_div(arb_mat_entry(X, i, c),
                    arb_mat_entry(X, i, c),
                    arb_mat_entry(L, i, i), prec);
        }

        /* solve Ux = y */
        for (i = n - 1; i >= 0; i--)
        {
            for (j = i + 1; j < n; j++)
            {
                arb_submul(arb_mat_entry(X, i, c),
                    arb_mat_entry(L, j, i), arb_mat_entry(X, j, c), prec);
            }
        }
    }
}
