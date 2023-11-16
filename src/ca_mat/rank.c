/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "ca_mat.h"

/* todo: different algorithms... */
int
ca_mat_rank(slong * rank, const ca_mat_t A, ca_ctx_t ctx)
{
    slong n, m;
    slong * P;
    int success;
    ca_mat_t T;

    n = ca_mat_nrows(A);
    m = ca_mat_ncols(A);

    if (n == 0 || m == 0)
    {
        *rank = 0;
        return 1;
    }

    ca_mat_init(T, n, m, ctx);
    P = _perm_init(n);
    success = ca_mat_lu(rank, P, T, A, 0, ctx);
    ca_mat_clear(T, ctx);
    _perm_clear(P);

    return success;
}
