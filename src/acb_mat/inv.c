/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_mat_inv(acb_mat_t X, const acb_mat_t A, slong prec)
{
    if (X == A)
    {
        int r;
        acb_mat_t T;
        acb_mat_init(T, acb_mat_nrows(A), acb_mat_ncols(A));
        r = acb_mat_inv(T, A, prec);
        acb_mat_swap(T, X);
        acb_mat_clear(T);
        return r;
    }

    acb_mat_one(X);
    return acb_mat_solve(X, A, X, prec);
}
