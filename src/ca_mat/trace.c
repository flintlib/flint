/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_trace(ca_t trace, const ca_mat_t mat, ca_ctx_t ctx)
{
    slong i;

    if (!ca_mat_is_square(mat))
    {
        flint_throw(FLINT_ERROR, "ca_mat_trace: a square matrix is required!\n");
    }

    if (ca_mat_is_empty(mat))
    {
        ca_zero(trace, ctx);
        return;
    }

    ca_set(trace, ca_mat_entry(mat, 0, 0), ctx);

    for (i = 1; i < ca_mat_nrows(mat); i++)
    {
        ca_add(trace, trace, ca_mat_entry(mat, i, i), ctx);
    }
}
