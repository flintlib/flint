/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "ca_mat.h"

void
ca_mat_set_fmpz_mat(ca_mat_t dest, const fmpz_mat_t src, ca_ctx_t ctx)
{
    slong i, j;

    if (ca_mat_ncols(dest) != 0)
    {
        for (i = 0; i < ca_mat_nrows(dest); i++)
            for (j = 0; j < ca_mat_ncols(dest); j++)
                ca_set_fmpz(ca_mat_entry(dest, i, j),
                    fmpz_mat_entry(src, i, j), ctx);
    }
}
