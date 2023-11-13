/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_mat.h"

void
d_mat_swap(d_mat_t mat1, d_mat_t mat2)
{
    FLINT_SWAP(d_mat_struct, *mat1, *mat2);
}

void
d_mat_swap_entrywise(d_mat_t mat1, d_mat_t mat2)
{
    slong i, j;
    for (i = 0; i < d_mat_nrows(mat1); i++)
    {
        double * row1 = mat1->rows[i];
        double * row2 = mat2->rows[i];
        for (j = 0; j < d_mat_ncols(mat1); j++)
            FLINT_SWAP(double, row1[j], row2[j]);
    }
}
