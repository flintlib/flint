/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_mat_solve_lu(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    int result;
    slong n, m, *perm;
    acb_mat_t LU;

    n = acb_mat_nrows(A);
    m = acb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    perm = _perm_init(n);
    acb_mat_init(LU, n, n);

    result = acb_mat_lu(perm, LU, A, prec);

    if (result)
        acb_mat_solve_lu_precomp(X, perm, LU, B, prec);

    acb_mat_clear(LU);
    _perm_clear(perm);

    return result;
}
