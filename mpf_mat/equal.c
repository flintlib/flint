/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpf_mat.h"

int
mpf_mat_equal(const mpf_mat_t mat1, const mpf_mat_t mat2)
{
    slong j;

    if (mat1->r != mat2->r || mat1->c != mat2->c)
    {
        return 0;
    }

    if (mat1->r == 0 || mat1->c == 0)
        return 1;

    for (j = 0; j < mat1->r; j++)
    {
        if (!_mpf_vec_equal(mat1->rows[j], mat2->rows[j], mat1->c))
        {
            return 0;
        }
    }

    return 1;
}
