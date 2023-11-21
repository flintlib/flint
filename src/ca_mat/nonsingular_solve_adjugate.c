/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

truth_t
ca_mat_nonsingular_solve_adjugate(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
{
    truth_t result;
    ca_t det;
    ca_mat_t T;

    ca_init(det, ctx);
    ca_mat_init(T, ca_mat_nrows(A), ca_mat_ncols(A), ctx);
    ca_mat_adjugate(T, det, A, ctx);

    result = ca_check_is_zero(det, ctx);
    if (result == T_TRUE)
        result = T_FALSE;
    else if (result == T_FALSE)
        result = T_TRUE;

    if (result == T_TRUE)
    {
        ca_mat_mul(X, T, B, ctx);
        ca_mat_div_ca(X, X, det, ctx);
    }

    ca_mat_clear(T, ctx);
    ca_clear(det, ctx);
    return result;
}
