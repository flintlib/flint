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
ca_mat_add(ca_mat_t res,
        const ca_mat_t mat1, const ca_mat_t mat2, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(mat1); i++)
        for (j = 0; j < ca_mat_ncols(mat1); j++)
            ca_add(ca_mat_entry(res, i, j),
                ca_mat_entry(mat1, i, j),
                ca_mat_entry(mat2, i, j), ctx);
}
