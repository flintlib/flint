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
arb_mat_spd_inv(arb_mat_t X, const arb_mat_t A, slong prec)
{
    slong n;
    arb_mat_t L;
    int result;

    if (!arb_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "arb_mat_spd_inv: a square matrix is required\n");
    }

    if (arb_mat_nrows(X) != arb_mat_nrows(A) ||
        arb_mat_ncols(X) != arb_mat_ncols(A))
    {
        flint_throw(FLINT_ERROR, "arb_mat_spd_inv: incompatible dimensions\n");
    }

    if (arb_mat_is_empty(A))
        return 1;

    n = arb_mat_nrows(A);

    if (n == 1)
    {
        if (arb_is_positive(arb_mat_entry(A, 0, 0)))
        {
            arb_inv(arb_mat_entry(X, 0, 0), arb_mat_entry(A, 0, 0), prec);
            return 1;
        }
        else
        {
            return 0;
        }
    }

    arb_mat_init(L, n, n);
    arb_mat_set(L, A);

    if (_arb_mat_cholesky_banachiewicz(L, prec))
    {
        arb_mat_inv_cho_precomp(X, L, prec);
        result = 1;
    }
    else
    {
        result = 0;
    }

    arb_mat_clear(L);
    return result;
}
