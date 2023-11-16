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
ca_mat_nonsingular_fflu(slong * P, ca_mat_t LU, ca_t den, const ca_mat_t A, ca_ctx_t ctx)
{
    if (ca_mat_is_empty(A))
    {
        ca_one(den, ctx);
        return T_TRUE;
    }
    else
    {
        int success;
        slong rank;
        success = ca_mat_fflu(&rank, P, LU, den, A, 1, ctx);

        if (success == 0)
            return T_UNKNOWN;

        if (rank == 0)
            return T_FALSE;

        return T_TRUE;
    }
}
