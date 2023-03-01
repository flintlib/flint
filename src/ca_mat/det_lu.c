/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

int
ca_mat_det_lu(ca_t res, const ca_mat_t A, ca_ctx_t ctx)
{
    truth_t invertible;
    slong * P;
    ca_mat_t T;
    slong i, n;

    n = ca_mat_nrows(A);
    P = _perm_init(n);
    ca_mat_init(T, n, n, ctx);
    invertible = ca_mat_nonsingular_lu(P, T, A, ctx);

    if (invertible == T_FALSE)
    {
        ca_zero(res, ctx);
    }
    else if (invertible == T_TRUE)
    {
        ca_one(res, ctx);
        for (i = 0; i < n; i++)
            ca_mul(res, res, ca_mat_entry(T, i, i), ctx);
        if (_perm_parity(P, n))
            ca_neg(res, res, ctx);
    }
    else
    {
        ca_unknown(res, ctx);
    }

    ca_mat_clear(T, ctx);
    _perm_clear(P);
    return invertible != T_UNKNOWN;
}
