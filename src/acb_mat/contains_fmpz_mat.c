/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_mat_contains_fmpz_mat(const acb_mat_t mat1, const fmpz_mat_t mat2)
{
    slong i, j;

    if ((acb_mat_nrows(mat1) != acb_mat_nrows(mat2)) ||
        (acb_mat_ncols(mat1) != acb_mat_ncols(mat2)))
        return 0;

    for (i = 0; i < acb_mat_nrows(mat1); i++)
        for (j = 0; j < acb_mat_ncols(mat1); j++)
            if (!acb_contains_fmpz(acb_mat_entry(mat1, i, j),
                fmpz_mat_entry(mat2, i, j)))
                return 0;

    return 1;
}
