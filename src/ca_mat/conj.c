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
ca_mat_conj(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j;

    if ((ca_mat_nrows(B) != ca_mat_nrows(A)) ||
        (ca_mat_ncols(B) != ca_mat_ncols(A)))
    {
        flint_throw(FLINT_ERROR, "ca_mat_conj: incompatible dimensions.\n");
    }

    for (i = 0; i < ca_mat_nrows(A); i++)
        for (j = 0; j < ca_mat_ncols(A); j++)
            ca_conj(ca_mat_entry(B, i, j), ca_mat_entry(A, i, j), ctx);
}
