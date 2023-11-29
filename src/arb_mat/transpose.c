/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_transpose(arb_mat_t B, const arb_mat_t A)
{
    slong i, j;

    if (arb_mat_nrows(B) != arb_mat_ncols(A) || arb_mat_ncols(B) != arb_mat_nrows(A))
    {
        flint_throw(FLINT_ERROR, "Exception (arb_mat_transpose). Incompatible dimensions.\n");
    }

    if (arb_mat_is_empty(A))
        return;

    if (A == B)  /* In-place, guaranteed to be square */
    {
        for (i = 0; i < arb_mat_nrows(A) - 1; i++)
        {
            for (j = i + 1; j < arb_mat_ncols(A); j++)
            {
                arb_swap(arb_mat_entry(A, i, j), arb_mat_entry(A, j, i));
            }
        }
    }
    else  /* Not aliased; general case */
    {
        for (i = 0; i < arb_mat_nrows(B); i++)
            for (j = 0; j < arb_mat_ncols(B); j++)
                arb_set(arb_mat_entry(B, i, j), arb_mat_entry(A, j, i));
    }
}
