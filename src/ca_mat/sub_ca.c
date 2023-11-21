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
ca_mat_sub_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx)
{
    slong i, j, n;

    if (a == y)
    {
        n = FLINT_MIN(ca_mat_nrows(a), ca_mat_ncols(a));
        for (i = 0; i < n; i++)
            ca_sub(ca_mat_entry(y, i, i), ca_mat_entry(y, i, i), x, ctx);
    }
    else
    {
        for (i = 0; i < ca_mat_nrows(a); i++)
        {
            for (j = 0; j < ca_mat_ncols(a); j++)
            {
                if (i == j)
                    ca_sub(ca_mat_entry(y, i, j), ca_mat_entry(a, i, j), x, ctx);
                else
                    ca_set(ca_mat_entry(y, i, j), ca_mat_entry(a, i, j), ctx);
            }
        }
    }
}
