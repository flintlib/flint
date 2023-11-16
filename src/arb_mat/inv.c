/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
arb_mat_inv(arb_mat_t X, const arb_mat_t A, slong prec)
{
    if (X == A)
    {
        int r;
        arb_mat_t T;
        arb_mat_init(T, arb_mat_nrows(A), arb_mat_ncols(A));
        r = arb_mat_inv(T, A, prec);
        arb_mat_swap(T, X);
        arb_mat_clear(T);
        return r;
    }

    arb_mat_one(X);
    return arb_mat_solve(X, A, X, prec);
}
