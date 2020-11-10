/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

truth_t
ca_mat_inv(ca_mat_t X, const ca_mat_t A, ca_ctx_t ctx)
{
    slong n;
    truth_t success;
    ca_mat_t T;

    n = ca_mat_nrows(A);

    if (n == 0)
        return T_TRUE;

    if (n <= 4)
    {
        ca_t det;
        ca_init(det, ctx);
        ca_mat_adjugate_cofactor(X, det, A, ctx);

        success = ca_check_is_zero(det, ctx);
        if (success == T_FALSE)
        {
            ca_mat_div_ca(X, X, det, ctx);
            success = T_TRUE;
        }
        else if (success == T_TRUE)
            success = T_FALSE;

        ca_clear(det, ctx);
        return success;
    }

    ca_mat_init(T, n, n, ctx);
    ca_mat_one(T, ctx);
    success = ca_mat_nonsingular_solve(X, A, T, ctx);
    ca_mat_clear(T, ctx);
    return success;
}
