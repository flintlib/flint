/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_transfer(ca_mat_t res, ca_ctx_t res_ctx, const ca_mat_t src, ca_ctx_t src_ctx)
{
    slong i, j;

    if (res_ctx == src_ctx)
    {
        ca_mat_set(res, src, res_ctx);
    }
    else
    {
        for (i = 0; i < ca_mat_nrows(src); i++)
            for (j = 0; j < ca_mat_ncols(src); j++)
                ca_transfer(ca_mat_entry(res, i, j), res_ctx, ca_mat_entry(src, i, j), src_ctx);
    }
}
