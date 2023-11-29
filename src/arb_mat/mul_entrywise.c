/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_mul_entrywise(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong i, j;

    if (arb_mat_nrows(A) != arb_mat_nrows(B) ||
        arb_mat_ncols(A) != arb_mat_ncols(B))
    {
        flint_throw(FLINT_ERROR, "arb_mat_mul_entrywise: incompatible dimensions\n");
    }

    for (i = 0; i < arb_mat_nrows(A); i++)
    {
        for (j = 0; j < arb_mat_ncols(A); j++)
        {
            arb_mul(arb_mat_entry(C, i, j),
                    arb_mat_entry(A, i, j),
                    arb_mat_entry(B, i, j), prec);
        }
    }
}
