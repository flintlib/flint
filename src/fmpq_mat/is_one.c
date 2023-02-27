/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

int
fmpq_mat_is_one(const fmpq_mat_t mat)
{
    slong i, j;

    if (mat->r == 0 || mat->c == 0)
        return 1;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            if (fmpq_cmp_ui(fmpq_mat_entry(mat, i, j), i == j) != 0)
            {
                return 0;
            }
        }
    }

    return 1;
}
