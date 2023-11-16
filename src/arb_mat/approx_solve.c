/*
    Copyright (C) 2018 arbguest

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "arb_mat.h"

int
arb_mat_approx_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong n, m, *perm;
    arb_mat_t LU;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    perm = _perm_init(n);
    arb_mat_init(LU, n, n);

    result = arb_mat_approx_lu(perm, LU, A, prec);

    if (result)
        arb_mat_approx_solve_lu_precomp(X, perm, LU, B, prec);

    arb_mat_clear(LU);
    _perm_clear(perm);

    return result;
}
