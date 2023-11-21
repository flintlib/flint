/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"
#include "arb_mat.h"

int
arb_mat_contains_fmpq_mat(const arb_mat_t mat1, const fmpq_mat_t mat2)
{
    slong i, j;

    if ((arb_mat_nrows(mat1) != arb_mat_nrows(mat2)) ||
        (arb_mat_ncols(mat1) != arb_mat_ncols(mat2)))
        return 0;

    for (i = 0; i < arb_mat_nrows(mat1); i++)
        for (j = 0; j < arb_mat_ncols(mat1); j++)
            if (!arb_contains_fmpq(arb_mat_entry(mat1, i, j),
                fmpq_mat_entry(mat2, i, j)))
                return 0;

    return 1;
}
