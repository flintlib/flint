/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_mul_entrywise(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong i, j;

    if (acb_mat_nrows(A) != acb_mat_nrows(B) ||
        acb_mat_ncols(A) != acb_mat_ncols(B))
    {
        flint_throw(FLINT_ERROR, "acb_mat_mul_entrywise: incompatible dimensions\n");
    }

    for (i = 0; i < acb_mat_nrows(A); i++)
    {
        for (j = 0; j < acb_mat_ncols(A); j++)
        {
            acb_mul(acb_mat_entry(C, i, j),
                    acb_mat_entry(A, i, j),
                    acb_mat_entry(B, i, j), prec);
        }
    }
}
