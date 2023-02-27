/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void fmpq_mat_one(fmpq_mat_t mat)
{
    slong i, j, min;

    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            fmpq_zero(fmpq_mat_entry(mat, i, j));

    min = FLINT_MIN(mat->r, mat->c);

    for (i = 0; i < min; i++)
        fmpq_one(fmpq_mat_entry(mat, i, i));
}

