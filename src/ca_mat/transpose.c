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
ca_mat_transpose(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j;

    if (ca_mat_nrows(B) != ca_mat_ncols(A) || ca_mat_ncols(B) != ca_mat_nrows(A))
    {
        flint_throw(FLINT_ERROR, "Exception (ca_mat_transpose). Incompatible dimensions.\n");
    }

    if (ca_mat_is_empty(A))
        return;

    if (A == B)  /* In-place, guaranteed to be square */
    {
        for (i = 0; i < ca_mat_nrows(A) - 1; i++)
        {
            for (j = i + 1; j < ca_mat_ncols(A); j++)
            {
                ca_swap(ca_mat_entry(A, i, j), ca_mat_entry(A, j, i), ctx);
            }
        }
    }
    else  /* Not aliased; general case */
    {
        for (i = 0; i < ca_mat_nrows(B); i++)
            for (j = 0; j < ca_mat_ncols(B); j++)
                ca_set(ca_mat_entry(B, i, j), ca_mat_entry(A, j, i), ctx);
    }
}
