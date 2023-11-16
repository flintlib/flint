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
ca_mat_check_is_one(const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j;
    truth_t res, eq;

    res = T_TRUE;

    for (i = 0; i < ca_mat_nrows(A); i++)
    {
        for (j = 0; j < ca_mat_ncols(A); j++)
        {
            if (i == j)
                eq = ca_check_is_one(ca_mat_entry(A, i, j), ctx);
            else
                eq = ca_check_is_zero(ca_mat_entry(A, i, j), ctx);

            if (eq == T_FALSE)
                return T_FALSE;
            if (eq == T_UNKNOWN)
                res = T_UNKNOWN;
        }
    }

    return res;
}
