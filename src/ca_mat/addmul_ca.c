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
ca_mat_addmul_ca(ca_mat_t y, const ca_mat_t a, const ca_t x, ca_ctx_t ctx)
{
    slong i, j;

    ca_t t;
    ca_init(t, ctx);

    for (i = 0; i < ca_mat_nrows(a); i++)
    {
        for (j = 0; j < ca_mat_ncols(a); j++)
        {
            ca_mul(t, ca_mat_entry(a, i, j), x, ctx);
            ca_add(ca_mat_entry(y, i, j), ca_mat_entry(y, i, j), t, ctx);
        }
    }

    ca_clear(t, ctx);
}
