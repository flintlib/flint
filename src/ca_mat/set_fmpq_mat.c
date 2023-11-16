/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"
#include "ca_mat.h"

void
ca_mat_set_fmpq_mat(ca_mat_t dest, const fmpq_mat_t src, ca_ctx_t ctx)
{
    slong i, j;

    if (ca_mat_ncols(dest) != 0)
    {
        for (i = 0; i < ca_mat_nrows(dest); i++)
            for (j = 0; j < ca_mat_ncols(dest); j++)
                ca_set_fmpq(ca_mat_entry(dest, i, j),
                    fmpq_mat_entry(src, i, j), ctx);
    }
}
