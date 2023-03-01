/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_neg(ca_mat_t dest, const ca_mat_t src, ca_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < ca_mat_nrows(src); i++)
        for (j = 0; j < ca_mat_ncols(src); j++)
            ca_neg(ca_mat_entry(dest, i, j),
                ca_mat_entry(src, i, j), ctx);
}
