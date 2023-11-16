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

truth_t
ca_mat_nonsingular_solve_fflu(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
{
    truth_t result;
    slong n, m, *perm;
    ca_mat_t LU;
    ca_t den;

    n = ca_mat_nrows(A);
    m = ca_mat_ncols(X);

    if (n == 0)
        return T_TRUE;

    perm = _perm_init(n);
    ca_mat_init(LU, n, n, ctx);
    ca_init(den, ctx);

    result = ca_mat_nonsingular_fflu(perm, LU, den, A, ctx);

    if (result == T_TRUE && m != 0)
        ca_mat_solve_fflu_precomp(X, perm, LU, den, B, ctx);

    ca_mat_clear(LU, ctx);
    _perm_clear(perm);
    ca_clear(den, ctx);

    return result;
}
