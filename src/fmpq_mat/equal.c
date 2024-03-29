/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

int
fmpq_mat_equal(const fmpq_mat_t mat1, const fmpq_mat_t mat2)
{
    slong i, j;

    if (mat1->r != mat2->r || mat1->c != mat2->c)
        return 0;

    for (i = 0; i < mat1->r; i++)
    {
        for (j = 0; j < mat1->c; j++)
        {
            if (!fmpq_equal(fmpq_mat_entry(mat1, i, j),
                            fmpq_mat_entry(mat2, i, j)))
                return 0;
        }
    }

    return 1;
}
