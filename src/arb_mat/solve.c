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
arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong n = arb_mat_nrows(A);

    if (n <= 4 || prec > 10.0 * n)
        return arb_mat_solve_lu(X, A, B, prec);
    else
        return arb_mat_solve_precond(X, A, B, prec);
}
