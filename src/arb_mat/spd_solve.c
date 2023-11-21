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
arb_mat_spd_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong n, m;
    arb_mat_t L;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    arb_mat_init(L, n, n);

    result = arb_mat_cho(L, A, prec);

    if (result)
        arb_mat_solve_cho_precomp(X, L, B, prec);

    arb_mat_clear(L);

    return result;
}
