/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
_arb_mat_cholesky_banachiewicz(arb_mat_t A, slong prec)
{
    slong n, i, j, k;

    n = arb_mat_nrows(A);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            for (k = 0; k < j; k++)
            {
                arb_submul(arb_mat_entry(A, i, j),
                           arb_mat_entry(A, i, k),
                           arb_mat_entry(A, j, k), prec);
            }
            arb_div(arb_mat_entry(A, i, j),
                    arb_mat_entry(A, i, j),
                    arb_mat_entry(A, j, j), prec);
        }
        for (k = 0; k < i; k++)
        {
            arb_submul(arb_mat_entry(A, i, i),
                       arb_mat_entry(A, i, k),
                       arb_mat_entry(A, i, k), prec);
        }
        if (!arb_is_positive(arb_mat_entry(A, i, i)))
            return 0;
        arb_sqrt(arb_mat_entry(A, i, i),
                 arb_mat_entry(A, i, i), prec);
    }

    return 1;
}

int
arb_mat_cho(arb_mat_t L, const arb_mat_t A, slong prec)
{
    slong n;

    if (!arb_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "arb_mat_cho: a square matrix is required\n");
    }

    if (arb_mat_nrows(L) != arb_mat_nrows(A) ||
        arb_mat_ncols(L) != arb_mat_ncols(A))
    {
        flint_throw(FLINT_ERROR, "arb_mat_cho: incompatible dimensions\n");
    }

    if (arb_mat_is_empty(A))
        return 1;

    n = arb_mat_nrows(A);

    if (n == 1)
    {
        if (arb_is_positive(arb_mat_entry(A, 0, 0)))
        {
            arb_sqrt(arb_mat_entry(L, 0, 0), arb_mat_entry(A, 0, 0), prec);
            return 1;
        }
        else
        {
            return 0;
        }
    }

    arb_mat_set(L, A);

    if (!_arb_mat_cholesky_banachiewicz(L, prec))
        return 0;

    /* set the strictly upper triangular region of L to zero */
    {
        slong i, j;
        for (i = 0; i < n; i++)
            for (j = i+1; j < n; j++)
                arb_zero(arb_mat_entry(L, i, j));
    }

    return 1;
}
