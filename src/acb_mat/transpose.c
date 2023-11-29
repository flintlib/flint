/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_transpose(acb_mat_t B, const acb_mat_t A)
{
    slong i, j;

    if (acb_mat_nrows(B) != acb_mat_ncols(A) || acb_mat_ncols(B) != acb_mat_nrows(A))
    {
        flint_throw(FLINT_ERROR, "Exception (acb_mat_transpose). Incompatible dimensions.\n");
    }

    if (acb_mat_is_empty(A))
        return;

    if (A == B)  /* In-place, guaranteed to be square */
    {
        for (i = 0; i < acb_mat_nrows(A) - 1; i++)
        {
            for (j = i + 1; j < acb_mat_ncols(A); j++)
            {
                acb_swap(acb_mat_entry(A, i, j), acb_mat_entry(A, j, i));
            }
        }
    }
    else  /* Not aliased; general case */
    {
        for (i = 0; i < acb_mat_nrows(B); i++)
            for (j = 0; j < acb_mat_ncols(B); j++)
                acb_set(acb_mat_entry(B, i, j), acb_mat_entry(A, j, i));
    }
}

