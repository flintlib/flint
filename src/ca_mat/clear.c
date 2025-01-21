/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_vec.h"
#include "ca_mat.h"

void
ca_mat_clear(ca_mat_t mat, ca_ctx_t ctx)
{
    if (mat->entries != NULL)
    {
        slong i, j;
        for (i = 0; i < mat->r; i++)
            for (j = 0; j < mat->c; j++)
                ca_clear(ca_mat_entry(mat, i, j), ctx);
        flint_free(mat->entries);
    }
}
