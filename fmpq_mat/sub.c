/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void fmpq_mat_sub(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            fmpq_sub(fmpq_mat_entry(mat, i, j), 
                     fmpq_mat_entry(mat1, i, j), fmpq_mat_entry(mat2, i, j));
}

